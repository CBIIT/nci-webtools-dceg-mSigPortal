import { Form, Row, Col, Button } from 'react-bootstrap';
import { useForm, Controller } from 'react-hook-form';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faFolderMinus } from '@fortawesome/free-solid-svg-icons';
import Select from '../../../controls/select/selectForm';
import { useSelector, useDispatch } from 'react-redux';
import { actions as exposureActions } from '../../../../services/store/exposure';
import { actions as modalActions } from '../../../../services/store/modal';
import {
  resetExplorationApi,
  useExposureOptionsQuery,
} from '../../../../services/store/rootApi';
import {
  useUploadExplorationMutation,
  useSubmitExplorationMutation,
} from './apiSlice';

const actions = { ...exposureActions, ...modalActions };
const { Group, Check, Label } = Form;

export default function PublicForm() {
  const dispatch = useDispatch();
  const mergeForm = async (state) =>
    await dispatch(actions.mergeExposure({ userForm: state }));
  const mergeMain = async (state) =>
    await dispatch(actions.mergeExposure({ main: state }));
  const resetExposure = (_) => dispatch(actions.resetExposure());
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const {
    submitted,
    loading,
    userNameOptions,
    userSampleOptions,
    id,
    usePublicSignature,
  } = useSelector((state) => state.exposure.main);

  const {
    data: exposureOptions,
    isFetching,
    isError,
  } = useExposureOptionsQuery();
  const [handleUpload, { isLoading: isUploading }] =
    useUploadExplorationMutation();
  const [submitExploration, { isLoading: loadingUserExposure }] =
    useSubmitExplorationMutation();

  const defaultFormValues = {
    exposureFile: '',
    matrixFile: '',
    signatureFile: '',
    study: { label: 'PCAWG', value: 'PCAWG' },
    signatureSetName: {
      label: 'COSMIC_v3_Signatures_GRCh37_SBS96',
      value: 'COSMIC_v3_Signatures_GRCh37_SBS96',
    },
    genome: { label: 'GRCh37', value: 'GRCh37' },
  };

  const {
    control,
    reset: resetForm,
    resetField,
    handleSubmit,
    watch,
    getValues,
    setValue,
    formState: { errors },
  } = useForm({ defaultValues: defaultFormValues });

  const formStudy = watch('study');
  const { exposureFile, matrixFile, signatureFile } = watch();

  // call calculate after receiving signature and sample name options
  // useEffect(() => {
  //   if (
  //     !submitted &&
  //     !loading &&
  //     userNameOptions.length &&
  //     userSampleOptions.length
  //   )
  //     calculate('all', id);
  // }, [userNameOptions, userSampleOptions, loading]);

  async function loadExample(type) {
    const filepath = `assets/exampleInput/Sherlock_SBS96_${type}.txt`;
    const filename = filepath.split('/').slice(-1)[0];
    if (`${type}File` != filename) {
      if (type == 'exposure') {
        setValue(
          'exposureFile',
          new File([await (await fetch(filepath)).blob()], filename)
        );
      } else if (type == 'matrix') {
        setValue(
          'matrixFile',
          new File([await (await fetch(filepath)).blob()], filename)
        );
      } else if (type == 'signature') {
        setValue(
          'signatureFile',
          new File([await (await fetch(filepath)).blob()], filename)
        );
      }
    }
  }

  //  get signature and sample names. useEffect will call main calculate function
  async function handleCalculate() {
    mergeForm({ loading: true });
    try {
      const { projectID: id, exposureData } = await handleUpload();
      // get signature name options, ignore sample key
      const nameOptions = Object.keys(exposureData[0]).filter(
        (key) => key != 'Samples'
      );
      const sampleOptions = [
        ...new Set(exposureData.map(({ Samples }) => Samples)),
      ];

      dispatch(
        actions.mergeExposure({
          main: {
            id,
            userNameOptions: nameOptions,
            userSampleOptions: sampleOptions,
          },
          msIndividual: { sample: sampleOptions[0] },
          msAssociation: {
            signatureName1: nameOptions[0],
            signatureName2: nameOptions[1],
          },
          msBurden: { signatureName: nameOptions[0] },
        })
      );
    } catch (err) {
      mergeError(err);
    }
    mergeForm({ loading: false });
  }

  const studyOptions = exposureOptions
    ? [...new Set(exposureOptions.map((e) => e.study))].sort().map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  const signatureSetOptions = (study) => {
    if (exposureOptions && study.value) {
      return [
        ...new Set(
          exposureOptions
            .filter((e) => e.study == study.value)
            .map((e) => e.signatureSetName)
        ),
      ]
        .sort()
        .map((e) => ({
          label: e,
          value: e,
        }));
    } else {
      return [];
    }
  };

  function handleStudy(study) {
    const signatureSets = signatureSetOptions(study);

    setValue('study', study);
    setValue('signatureSetName', signatureSets[0]);
  }

  async function onSubmit(data) {
    try {
      mergeMain({ submitted: true, loading: true });
      mergeForm(data);

      // create form data for upload
      const formData = new FormData();
      formData.append('exposureFile', exposureFile);
      formData.append('matrixFile', matrixFile);
      if (!usePublicSignature) formData.append('signatureFile', signatureFile);
      const { projectID: id } = await handleUpload(formData).unwrap();

      // submit after upload
      await submitExploration({
        id,
        exposureFile: data.exposureFile.name,
        matrixFile: data.matrixFile.name,
        signatureFile: data?.signatureFile.name,
      }).unwrap();

      // // get signature name options, ignore sample key
      // const nameOptions = Object.keys(exposureData[0]).filter(
      //   (key) => key != 'Samples'
      // );
      // const sampleOptions = [
      //   ...new Set(exposureData.map(({ Samples }) => Samples)),
      // ];

      // dispatch(
      //   actions.mergeExposure({
      //     main: {
      //       id,
      //       userNameOptions: nameOptions,
      //       userSampleOptions: sampleOptions,
      //     },
      //     msIndividual: { sample: sampleOptions[0] },
      //     msAssociation: {
      //       signatureName1: nameOptions[0],
      //       signatureName2: nameOptions[1],
      //     },
      //     msBurden: { signatureName: nameOptions[0] },
      //   })
      // );

      mergeMain({
        displayTab: 'tmb',
        id,
        samples: [1],
        signatureNames: [1],
      });
    } catch (error) {
      mergeError(error.message || 'Failed to submit');
    } finally {
      mergeMain({ loading: false });
    }
  }

  function handleReset() {
    window.location.hash = '#/exploration';
    resetForm();
    resetExposure();
    dispatch(resetExplorationApi);
    mergeMain({ source: 'user' });
  }

  return (
    <Form onSubmit={handleSubmit(onSubmit)}>
      <Row>
        <Col>
          <Group>
            <Label>Upload Exposure File</Label>
            <Row className="m-0">
              <Col lg="6" className="p-0">
                <Button
                  className="p-0 font-14"
                  disabled={submitted}
                  variant="link"
                  href={'assets/exampleInput/Sherlock_SBS96_exposure.txt'}
                  download
                >
                  Download Example
                </Button>
              </Col>
              <Col lg="6" className="p-0 d-flex">
                <Button
                  className={`p-0 ml-auto font-14`}
                  disabled={submitted || loading.active}
                  variant="link"
                  type="button"
                  onClick={() => loadExample('exposure')}
                >
                  Load Example
                </Button>
              </Col>
            </Row>
            <div className="d-flex">
              <Controller
                name="exposureFile"
                control={control}
                rules={{ required: true }}
                render={({ field }) => (
                  <Form.File
                    {...field}
                    value={''}
                    disabled={loading || submitted}
                    id="exposureFile"
                    title={exposureFile.name || 'Upload Exposure File'}
                    label={exposureFile.name || 'Upload Exposure File'}
                    // accept=".txt"
                    isInvalid={errors.exposureFile}
                    feedback="Upload an exposure file"
                    onChange={(e) => {
                      if (e.target.files.length) {
                        setValue('exposureFile', e.target.files[0]);
                      }
                    }}
                    custom
                  />
                )}
              />
              {exposureFile.size > 0 && (
                <Button
                  className="ml-1"
                  size="sm"
                  title="remove exposure file"
                  variant="danger"
                  disabled={loading || submitted}
                  onClick={() => resetField('exposureFile')}
                >
                  <FontAwesomeIcon icon={faFolderMinus} size="lg" />
                </Button>
              )}
            </div>
          </Group>
        </Col>
      </Row>
      <Row>
        <Col>
          <Group>
            <Label>Upload Matrix File</Label>
            <Row className="m-0">
              <Col lg="6" className="p-0">
                <Button
                  className="p-0 font-14"
                  disabled={submitted}
                  variant="link"
                  href={'assets/exampleInput/Sherlock_SBS96_matrix.txt'}
                  download
                >
                  Download Example
                </Button>
              </Col>
              <Col lg="6" className="p-0 d-flex">
                <Button
                  className={`p-0 ml-auto font-14`}
                  disabled={submitted || loading}
                  variant="link"
                  type="button"
                  onClick={() => loadExample('matrix')}
                >
                  Load Example
                </Button>
              </Col>
            </Row>
            <div className="d-flex">
              <Controller
                name="matrixFile"
                control={control}
                rules={{ required: true }}
                render={({ field }) => (
                  <Form.File
                    {...field}
                    value={''}
                    disabled={loading || submitted}
                    id="uploadMatrix"
                    title={matrixFile.name || 'Upload Matrix File'}
                    label={matrixFile.name || 'Upload Matrix File'}
                    // accept=".txt"
                    isInvalid={errors.matrixFile}
                    feedback="Upload a matrix file"
                    onChange={(e) => {
                      if (e.target.files.length) {
                        setValue('matrixFile', e.target.files[0]);
                      }
                    }}
                    custom
                  />
                )}
              />
              {matrixFile.size > 0 && (
                <Button
                  className="ml-1"
                  size="sm"
                  title="remove matrix file"
                  variant="danger"
                  disabled={loading || submitted}
                  onClick={() => resetField('matrixFile')}
                >
                  <FontAwesomeIcon icon={faFolderMinus} size="lg" />
                </Button>
              )}
            </div>
          </Group>
        </Col>
      </Row>
      <Row>
        <Col>
          <Group controlId="toggleSignatureSource" className="d-flex">
            <Label className="mr-4">Use Public Signature Data</Label>
            <Check inline id="toggleSignatureSource">
              <Check.Input
                disabled={loading || submitted}
                type="checkbox"
                value={usePublicSignature}
                checked={usePublicSignature}
                onChange={() =>
                  mergeMain({
                    usePublicSignature: !usePublicSignature,
                  })
                }
              />
            </Check>
          </Group>
        </Col>
      </Row>

      {usePublicSignature ? (
        <div>
          <Select
            name="study"
            label="Study"
            control={control}
            disabled={loading || submitted || isFetching}
            options={studyOptions}
            onChange={handleStudy}
          />
          <Select
            name="signatureSetName"
            label="Reference Signature Set"
            control={control}
            disabled={loading || submitted || isFetching}
            options={signatureSetOptions(formStudy)}
          />
        </div>
      ) : (
        <Row>
          <Col>
            <Group>
              <Label>Upload Signature Data</Label>
              <Row className="m-0">
                <Col lg="6" className="p-0">
                  <Button
                    className="p-0 font-14"
                    disabled={submitted}
                    variant="link"
                    href={'assets/exampleInput/Sherlock_SBS96_signature.txt'}
                    download
                  >
                    Download Example
                  </Button>
                </Col>
                <Col lg="6" className="p-0 d-flex">
                  <Button
                    className={`p-0 ml-auto font-14`}
                    disabled={submitted || loading.active}
                    variant="link"
                    type="button"
                    onClick={() => loadExample('signature')}
                  >
                    Load Example
                  </Button>
                </Col>
              </Row>
              <div className="d-flex">
                <Controller
                  name="signatureFile"
                  control={control}
                  rules={{ required: usePublicSignature }}
                  render={({ field }) => (
                    <Form.File
                      {...field}
                      value={''}
                      disabled={loading || submitted}
                      id="uploadSignature"
                      title={signatureFile.name || 'Upload Signature File'}
                      label={signatureFile.name || 'Upload Signature File'}
                      // accept=".txt"
                      isInvalid={errors.signatureFile}
                      feedback="Upload a signature file"
                      onChange={(e) => {
                        if (e.target.files.length) {
                          setValue('signatureFile', e.target.files[0]);
                        }
                      }}
                      custom
                    />
                  )}
                />
                {signatureFile.size > 0 && (
                  <Button
                    className="ml-1"
                    size="sm"
                    title="remove signature file"
                    variant="danger"
                    disabled={loading || submitted}
                    onClick={() => resetField('signatureFile')}
                  >
                    <FontAwesomeIcon icon={faFolderMinus} size="lg" />
                  </Button>
                )}
              </div>
            </Group>
          </Col>
        </Row>
      )}
      <Row>
        <Col>
          <Group>
            <Select
              name="genome"
              label="Genome"
              control={control}
              disabled={id}
              options={[
                { label: 'GRCh37', value: 'GRCh37' },
                { label: 'GRCh38', value: 'GRCh38' },
                { label: 'mm10', value: 'mm10' },
              ]}
            />
          </Group>
        </Col>
      </Row>
      <Row>
        <Col lg="6">
          <Button
            className="w-100 mb-3"
            variant="secondary"
            onClick={() => handleReset()}
          >
            Reset
          </Button>
        </Col>
        <Col lg="6">
          <Button
            disabled={
              loading || isFetching || isUploading || loadingUserExposure || id
            }
            className="w-100"
            variant="primary"
            type="submit"
          >
            Calculate
          </Button>
        </Col>
      </Row>
    </Form>
  );
}
