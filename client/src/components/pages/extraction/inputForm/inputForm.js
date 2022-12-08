import { useState } from 'react';
import { Form, Button, Row, Col } from 'react-bootstrap';
import { useForm, Controller } from 'react-hook-form';
import Select from '../../../controls/select/selectForm';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { useSelector, useDispatch } from 'react-redux';
import { actions as extractionActions } from '../../../../services/store/extraction';
import { actions as modalActions } from '../../../../services/store/modal';
import {
  resetExtractionApi,
  useSeqmatrixOptionsQuery,
} from '../../../../services/store/rootApi';
import { useSignatureOptionsQuery } from './apiSlice';

const actions = { ...extractionActions, ...modalActions };

export default function InputForm() {
  const store = useSelector((state) => state.extraction);
  const { submitted } = store.main;
  const { inputFilename } = store.inputForm;

  const dispatch = useDispatch();
  const mergeState = (state) =>
    dispatch(actions.mergeExtraction({ publicForm: state }));
  const resetExtraction = (_) => dispatch(actions.resetExtraction());
  const mergeMain = (state) =>
    dispatch(actions.mergeExtraction({ main: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  // query options to populate form
  const {
    data: seqmatrixOptions,
    error: seqmatrixError,
    isFetching: fetchingSeqmatrixOptions,
  } = useSeqmatrixOptionsQuery();
  const {
    data: signatureOptions,
    error: signatureError,
    isFetching: fetchingSignatureOptions,
  } = useSignatureOptionsQuery();

  // toggle visibility of advanced menu
  const [showAdvanced, setShowAdvanced] = useState(false);

  const genomeOptions = [
    { label: 'GRCh37', value: 'GRCh37' },
    { label: 'GRCh38', value: 'GRCh38' },
    { label: 'mm10', value: 'mm10' },
  ];
  const defaultValues = {
    source: 'public',
    study: { label: 'PCAWG', value: 'PCAWG' },
    cancer: { label: 'Lung-AdenoCA', value: 'Lung-AdenoCA' },
    strategy: { label: 'WGS', value: 'WGS' },
    genome: genomeOptions[0],
    exome: false,
    signatureSetName: {
      label: 'COSMIC_v3.3_Signatures_GRCh37_SBS96',
      value: 'COSMIC_v3.3_Signatures_GRCh37_SBS96',
    },
    signatureName: { label: 'all', value: 'all' },
    extractTool: {
      label: 'SigProfilerExtractor',
      value: 'SigProfilerExtractor',
    },
  };

  const {
    control,
    register,
    handleSubmit,
    reset: resetForm,
    setValue,
    watch,
    formState: { errors },
  } = useForm({ defaultValues });

  const {
    source,
    study,
    cancer,
    inputFile,
    exome,
    signatureSetName,
    useQueue,
  } = watch();

  function handleReset() {
    window.location.hash = '#/extraction';
    resetForm();
    resetExtraction();
    dispatch(resetExtractionApi);
  }

  async function onSubmit(data) {
    console.log(data);
    return;
    try {
      mergeMain({ submitted: true, loading: { active: true } });
      mergeState(data);
      const params = {
        study: data.study.value,
        cancer: data.cancer.value,
        strategy: data.strategy.value,
      };

      // let matrixData = [];
      // for await (const data of paginateQuery(fetchMatrix, params)) {
      //   matrixData = [...matrixData, ...data];
      // }
      const matrixData = []; //await fetchMatrix(params).unwrap();

      mergeMain({ matrixData, projectID: crypto.randomUUID() });
    } catch (error) {
      console.log(error);
      if (error.originalStatus == 504) {
        mergeMain({
          error: 'Please Reset Your Parameters and Try again.',
        });
        mergeError({
          visible: true,
          message:
            'Your submission has timed out. Please try again by submitting this job to a queue instead.',
        });
      } else {
        mergeMain({
          error: 'Please Reset Your Parameters and Try again.',
        });
        mergeError(error.data);
      }
    }
    mergeMain({ loading: { active: false } });
  }

  const studyOptions = seqmatrixOptions
    ? [...new Set(seqmatrixOptions.map((e) => e.study))].sort().map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  const cancerOptions = (study) => {
    if (seqmatrixOptions && study.value) {
      const options = [
        ...[
          ...new Set(
            seqmatrixOptions
              .filter((e) => e.study == study.value)
              .map((e) => e.cancer)
          ),
        ]
          .sort()
          .map((e) => ({
            label: e,
            value: e,
          })),
      ];
      if (study.value == 'PCAWG' || study.value == 'TCGA') {
        return [
          ...options,
          {
            label: 'PanCancer',
            value: '*ALL',
          },
        ];
      } else {
        return options;
      }
    } else {
      return [];
    }
  };

  const strategyOptions = (study, cancer) => {
    if (seqmatrixOptions && study.value) {
      const options =
        cancer.value == '*ALL'
          ? [
              ...new Set(
                seqmatrixOptions
                  .filter((e) => e.study == study.value)
                  .map((e) => e.strategy)
              ),
            ]
          : [
              ...new Set(
                seqmatrixOptions
                  .filter(
                    (e) => e.study == study.value && e.cancer == cancer.value
                  )
                  .map((e) => e.strategy)
              ),
            ];
      return options.sort().map((e) => ({
        label: e,
        value: e,
      }));
    } else {
      return [];
    }
  };

  function handleStudyChange(study) {
    const cancers = cancerOptions(study);
    const strategies = strategyOptions(study, cancers[0]);

    setValue('study', study);
    setValue('cancer', cancers[0]);
    setValue('strategy', strategies[0]);
  }

  function handleCancerChange(cancer) {
    const strategies = strategyOptions(study, cancer);

    setValue('cancer', cancer);
    setValue('strategy', strategies[0]);
  }

  const signatureSetOptions = signatureOptions
    ? [...new Set(signatureOptions.map((e) => e.signatureSetName))].map(
        (e) => ({ label: e, value: e })
      )
    : [];

  const signatureNameOptions = (signatureSetName) =>
    signatureOptions
      ? [
          { label: 'all', value: 'all' },
          ...[
            ...new Set(
              signatureOptions
                .filter((e) => e.signatureSetName == signatureSetName.value)
                .map((e) => e.signatureName)
            ),
          ].map((e) => ({ label: e, value: e })),
        ]
      : [];

  return (
    <Form onSubmit={handleSubmit(onSubmit)}>
      <LoadingOverlay active={fetchingSeqmatrixOptions} />
      {(seqmatrixError || signatureError) && (
        <p>There was an error retrieving public data options</p>
      )}
      <div className="border rounded p-2 mb-3">
        <Form.Group>
          <Form.Label className="mr-4">Data Source</Form.Label>
          <Controller
            name="source"
            control={control}
            render={({ field }) => (
              <Form.Check
                {...field}
                inline
                id="radioPublic"
                type="radio"
                label={<span className="font-weight-normal">Public</span>}
                value={'public'}
                checked={source == 'public'}
                disabled={submitted}
              />
            )}
          />
          <Controller
            name="source"
            control={control}
            render={({ field }) => (
              <Form.Check
                {...field}
                inline
                id="radioUser"
                type="radio"
                label={<span className="font-weight-normal">User</span>}
                value={'user'}
                checked={source == 'user'}
                disabled={submitted}
              />
            )}
          />
        </Form.Group>
        {source == 'public' ? (
          <div>
            <Select
              className="mb-2"
              name="study"
              label="Study"
              disabled={submitted || fetchingSeqmatrixOptions}
              options={studyOptions}
              control={control}
              onChange={handleStudyChange}
            />
            <Select
              className="mb-2"
              name="cancer"
              label="Cancer Type or Group"
              disabled={submitted || fetchingSeqmatrixOptions}
              options={cancerOptions(study)}
              control={control}
              onChange={handleCancerChange}
            />
            <Select
              className="mb-2"
              name="strategy"
              label="Experimental Strategy"
              disabled={submitted || fetchingSeqmatrixOptions}
              options={strategyOptions(study, cancer)}
              control={control}
            />
          </div>
        ) : (
          <div>
            <Select
              name="dataType"
              label="Data Type"
              disabled={submitted}
              options={[{ label: 'default', value: 'default' }]}
              control={control}
            />
            <Form.Group>
              <Form.Label>
                Upload File <span style={{ color: 'crimson' }}>*</span>
              </Form.Label>
              <Controller
                name="inputFile"
                control={control}
                rules={{ required: source == 'user' }}
                render={({ field }) => (
                  <Form.File
                    {...field}
                    value={''} // set dummy value for file input
                    disabled={submitted}
                    id="inputFile"
                    title={inputFile?.name || 'Upload Data File...'}
                    label={
                      inputFile?.name || inputFilename || 'Upload Data File...'
                    }
                    accept=".csv, .tsv, .vcf, .gz, .tar, .tar.gz, .txt"
                    isInvalid={errors.inputFile}
                    feedback="Please upload a data file"
                    onChange={(e) => {
                      if (e.target.files.length) {
                        setValue('inputFile', e.target.files[0]);
                      }
                    }}
                    custom
                  />
                )}
              />
            </Form.Group>
          </div>
        )}
      </div>

      <Select
        name="genome"
        label="Referencee Genome Build"
        disabled={submitted}
        options={genomeOptions}
        control={control}
      />
      <Form.Group>
        <Controller
          name="exome"
          control={control}
          render={({ field }) => (
            <Form.Check
              {...field}
              id="exome"
              type="checkbox"
              label={'Exome'}
              checked={exome}
              disabled={submitted}
            />
          )}
        />
      </Form.Group>
      <Form.Group>
        <Form.Label>Context Type</Form.Label>
        <Form.Control {...register('contextType')} defaultValue="default" />
      </Form.Group>
      <Select
        name="signatureSetName"
        label="Referencee Signature Set"
        disabled={submitted}
        options={signatureSetOptions}
        control={control}
      />
      <Select
        name="signatureName"
        label="Included Signature Names"
        disabled={submitted}
        options={signatureNameOptions(signatureSetName)}
        control={control}
        isMulti
      />
      <Select
        name="extractTool"
        label="Extract Tool"
        disabled={submitted}
        options={[
          { label: 'SigProfilerExtractor', value: 'SigProfilerExtractor' },
        ]}
        control={control}
      />

      <Button
        className="p-0"
        variant="link"
        onClick={(_) => setShowAdvanced(!showAdvanced)}
      >
        Advanced Parameters +
      </Button>

      <div className={showAdvanced ? 'd-block' : 'd-none'}>
        <Form.Label>Execution</Form.Label>
        <Form.Group>
          <Form.Label>CPU</Form.Label>
          <Form.Control {...register('cpu')} type="number" defaultValue={2} />
        </Form.Group>
        <Form.Group>
          <Form.Check
            {...register('gpu')}
            id="gpu"
            label="GPU"
            defaultChecked={false}
          />
        </Form.Group>
        <Form.Group>
          <Form.Label>Batch Size</Form.Label>
          <Form.Control
            {...register('batchSize')}
            type="number"
            defaultValue={1}
          />
        </Form.Group>
        <hr />
        <Form.Label>NMF Replicates</Form.Label>
        <Form.Group>
          <Form.Label>Minimum Signatures</Form.Label>
          <Form.Control
            {...register('minSignatures')}
            type="number"
            defaultValue={1}
          />
        </Form.Group>
        <Form.Group>
          <Form.Label>Maximum Signatures</Form.Label>
          <Form.Control
            {...register('maxSignatures')}
            type="number"
            defaultValue={15}
          />
        </Form.Group>
        <Form.Group>
          <Form.Group>
            <Form.Label>NMF Replicates Size</Form.Label>
            <Form.Control
              {...register('minSignatures')}
              type="number"
              defaultValue={100}
            />
          </Form.Group>
          <Form.Check
            {...register('resample')}
            label="Resamples"
            id="resample"
            defaultChecked={true}
          />
        </Form.Group>
        <Form.Group>
          <Form.Label>Seed</Form.Label>
          <Form.Control {...register('seed')} />
        </Form.Group>
      </div>

      <hr className="mb-3" />
      <div>
        <Form.Group controlId="toggleQueue" className="d-flex">
          <Form.Label className="mr-auto">
            Submit this job to a Queue
          </Form.Label>{' '}
          <Form.Check inline>
            <Controller
              name="useQueue"
              control={control}
              render={({ field }) => (
                <Form.Check.Input
                  {...field}
                  type="checkbox"
                  // disabled={submitted}
                  disabled={true}
                />
              )}
            />
          </Form.Check>
        </Form.Group>
        <Form.Group controlId="email">
          <Controller
            name="email"
            control={control}
            rules={{
              required: useQueue,
              pattern:
                /^[a-zA-Z0-9.!#$%&'*+/=?^_`{|}~-]+@[a-zA-Z0-9-]+(?:\.[a-zA-Z0-9-]+)*$/,
            }}
            render={({ field }) => (
              <Form.Control
                {...field}
                aria-label="Enter Email"
                placeholder="Enter Email"
                size="sm"
                type="email"
                disabled={!useQueue || submitted}
                isInvalid={errors.email}
              />
            )}
          />
          <Form.Control.Feedback type="invalid">
            {errors.email && 'Please provide a valid email'}
          </Form.Control.Feedback>
          <Form.Text className="text-muted">
            <i>
              Note: If sending to queue, when computation is completed, a
              notification will be sent to the e-mail entered above.
            </i>
          </Form.Text>
        </Form.Group>
      </div>

      <Row>
        <Col>
          <Button
            disabled={fetchingSeqmatrixOptions}
            className="w-100"
            variant="secondary"
            onClick={() => handleReset()}
          >
            Reset
          </Button>
        </Col>
        <Col>
          <Button
            disabled={fetchingSeqmatrixOptions || submitted}
            className="w-100"
            variant="primary"
            type="submit"
          >
            Submit
          </Button>
        </Col>
      </Row>
    </Form>
  );
}
