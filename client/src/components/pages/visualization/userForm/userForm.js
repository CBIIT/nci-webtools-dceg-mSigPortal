import { useEffect } from 'react';
import {
  Form,
  Button,
  Row,
  Col,
  Popover,
  OverlayTrigger,
} from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faInfoCircle, faFolderMinus } from '@fortawesome/free-solid-svg-icons';
import { useForm, Controller } from 'react-hook-form';
import Select from '../../../controls/select/selectForm';
import { actions as visualizationActions } from '../../../../services/store/visualization';
import { actions as modalActions } from '../../../../services/store/modal';
import { resetVisualizationApi } from '../../../../services/store/rootApi';
import {
  useVisualizationUserUploadMutation,
  useSubmitQueueMutation,
  useProfilerExtractionMutation,
  useUserMatrixMutation,
} from './apiSlice';

const actions = { ...visualizationActions, ...modalActions };

const { Title, Content } = Popover;

export default function UserForm() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.visualization);
  const mergeState = (state) =>
    dispatch(actions.mergeVisualization({ userForm: state }));
  const mergeMain = (state) =>
    dispatch(actions.mergeVisualization({ main: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));
  const mergeSuccess = (msg) =>
    dispatch(actions.mergeModal({ success: { visible: true, message: msg } }));
  const resetVisualization = (_) => dispatch(actions.resetVisualization());

  const { inputFilename, bedFilename, projectID, ...userForm } = store.userForm;
  const { submitted, source } = store.main;

  const [handleUpload, { isLoading: isUploading }] =
    useVisualizationUserUploadMutation();
  const [handleSubmitQueue, { isLoading: isQueueing }] =
    useSubmitQueueMutation();
  const [profilerExtraction, { isLoading }] = useProfilerExtractionMutation();
  const [fetchMatrix, { isLoading: loadingMatrix }] = useUserMatrixMutation();

  const formatOptions = [
    { label: 'VCF', value: 'vcf', example: 'demo_input_multi.vcf.gz' },
    { label: 'MAF', value: 'maf', example: 'demo_input_multi.MAF.gz' },
    { label: 'CSV', value: 'csv', example: 'demo_input_multi.csv' },
    { label: 'TSV', value: 'tsv', example: 'demo_input_multi.tsv' },
    {
      label: 'CATALOG TSV',
      value: 'catalog_tsv',
      example: 'demo_input_catalog.tsv',
    },
    {
      label: 'CATALOG CSV',
      value: 'catalog_csv',
      example: 'demo_input_catalog.csv',
    },
  ];

  const genomeOptions = [
    { label: 'GRCh37', value: 'GRCh37' },
    { label: 'GRCh38', value: 'GRCh38' },
    { label: 'mm10', value: 'mm10' },
  ];

  const defaultFormValues = {
    inputFormat: formatOptions[0],
    inputFile: '',
    genome: genomeOptions[0],
    strategy: 'WGS',
    mutationSplit: false,
    filter: '',
    bedFile: '',
    collapse: false,
    useQueue: false,
    email: '',
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

  const inputFormat = watch('inputFormat');
  const inputFile = watch('inputFile');
  const bedFile = watch('bedFile');
  const strategy = watch('strategy');
  const mutationSplit = watch('mutationSplit');
  const filter = watch('filter');
  const useQueue = watch('useQueue');

  function handleReset() {
    window.location.hash = '#/visualization';
    resetForm(defaultFormValues);
    resetVisualization();
    dispatch(resetVisualizationApi);
    mergeMain({ source: 'user' });
  }

  function removeFile() {
    resetField('inputFile');
  }

  function removeBedFile() {
    resetField('bedFile');
  }

  function handleSelectFormat(format) {
    removeFile();
    removeBedFile();
    setValue('inputFormat', format);
  }

  async function loadExample() {
    const filename = getValues('inputFormat').example;
    if (inputFilename != filename) {
      const file = 'assets/exampleInput/' + filename;
      setValue(
        'inputFile',
        new File([await (await fetch(file)).blob()], filename)
      );
    }
  }

  async function loadBed() {
    const filename = 'demo_input_bed.bed';
    if (bedFilename != filename) {
      const file = 'assets/exampleInput/' + filename;
      setValue(
        'bedFile',
        new File([await (await fetch(file)).blob()], filename)
      );
    }
  }

  // upload files and continue with web or queue submit
  async function onSubmit(data) {
    try {
      mergeMain({ submitted: true });
      const formData = new FormData();
      formData.append('inputFile', data.inputFile);
      if (bedFile.size) formData.append('bedFile', data.bedFile);

      const { projectID, ...filePaths } = await handleUpload(formData).unwrap();

      mergeState({
        ...data,
        ...filePaths,
        projectID,
        inputFilename: data.inputFile.name,
        bedFilename: data.bedFile.name,
      });
    } catch (error) {
      mergeError('An error occured while uploading files. Please try again.');
    }
  }

  async function submitWeb(args) {
    try {
      mergeMain({ loading: { active: true } });
      const { stdout, stderr, ...rest } = await profilerExtraction(
        args
      ).unwrap();
      const matrixData = await fetchMatrix({ userId: projectID }).unwrap();

      mergeMain({ projectID, matrixData, ...rest });
    } catch (error) {
      if (error.originalStatus == 504) {
        mergeMain({
          error: 'Please Reset Your Parameters and Try again.',
        });
        mergeError(
          'Your submission has timed out. Please try again by submitting this job to a queue instead.'
        );
      } else {
        mergeMain({
          error: 'Please Reset Your Parameters and Try again.',
        });

        const message = Object.values(error.data);
        mergeError(message);
      }
    }
    mergeMain({ loading: { active: false } });
  }

  async function submitQueue(args) {
    try {
      await handleSubmitQueue({ args, state: { ...store } });

      mergeSuccess(
        `Your job was successfully submitted to the queue. You will recieve an email at ${userForm.email} with your results.`
      );
    } catch (err) {
      mergeError('Failed to submit to queue. Please Try Again.');
    }
  }

  // perform web or queue submit after file upload and form is merged
  useEffect(() => {
    if (source == 'user' && projectID) {
      const pythonArgs = {
        inputFormat: ['-f', userForm.inputFormat.value],
        inputFile: ['-i', userForm.filePath],
        projectID: ['-p', projectID],
        genomeAssemblyVersion: ['-g', userForm.genome.value],
        strategy: ['-t', userForm.strategy],
        outputDir: ['-o', projectID],
      };

      if (useQueue) submitQueue(pythonArgs);
      else submitWeb(pythonArgs);
    }
  }, [projectID]);

  const msPopover = (
    <Popover id="popover-basic">
      <Title as="h3">Mutation Split</Title>
      <Content>
        <p>
          For each sample, split mutations into different groups according to
          the “Filter” column in VCF/CSV/TSV file. Splitting operation uses the
          “;” as separator.
        </p>
      </Content>
    </Popover>
  );

  const csPopover = (
    <Popover id="popover-basic">
      <Title as="h3">Mutation Split</Title>
      <Content>
        <p>
          A new sample called “All_Sample” will be added to the result, which
          combines the mutations from all samples.
        </p>
      </Content>
    </Popover>
  );

  return (
    <Form onSubmit={handleSubmit(onSubmit)}>
      <Select
        name="inputFormat"
        label="Select File Format"
        control={control}
        disabled={submitted}
        options={formatOptions}
        onChange={handleSelectFormat}
      />
      <Form.Group>
        <Form.Label>
          Upload File <span style={{ color: 'crimson' }}>*</span>
        </Form.Label>
        <Row className="m-0">
          <Col lg="6" className="p-0">
            <Button
              className="p-0 font-14"
              disabled={submitted}
              variant="link"
              href={'assets/exampleInput/all_demo_inputs.zip'}
              download
            >
              Download Example Data
            </Button>
          </Col>
          <Col lg="6" className="p-0 d-flex">
            <Button
              className={`p-0 ml-auto font-14`}
              disabled={submitted || isUploading || isQueueing || isLoading}
              variant="link"
              type="button"
              onClick={() => loadExample()}
            >
              Load Example Data
            </Button>
          </Col>
        </Row>
        <Row>
          <Col>
            <div className="d-flex">
              <Controller
                name="inputFile"
                control={control}
                rules={{ required: true }}
                render={({ field }) => (
                  <Form.File
                    {...field}
                    value={''} // set dummy value for file input
                    disabled={
                      submitted || isUploading || isQueueing || isLoading
                    }
                    id="inputFile"
                    title={inputFile.name || 'Upload Data File...'}
                    label={
                      inputFile.name || inputFilename || 'Upload Data File...'
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
              {inputFile.size > 0 && (
                <Button
                  className="ml-1"
                  size="sm"
                  title="Remove"
                  variant="danger"
                  disabled={submitted || isUploading || isQueueing || isLoading}
                  onClick={removeFile}
                >
                  <FontAwesomeIcon icon={faFolderMinus} size="lg" />
                </Button>
              )}
            </div>
          </Col>
        </Row>
      </Form.Group>
      <Select
        name="genome"
        label="Select Reference Genome Build"
        control={control}
        disabled={
          submitted || ['catalog_csv', 'catalog_tsv'].includes(inputFormat)
        }
        options={genomeOptions}
      />
      <Form.Group className="d-flex">
        <Form.Label className="mr-auto">Experimental Strategy</Form.Label>
        <Form.Check inline id="radioWGS">
          <Controller
            name="strategy"
            control={control}
            render={({ field }) => (
              <Form.Check.Input
                {...field}
                type="radio"
                value="WGS"
                checked={strategy == 'WGS'}
                disabled={
                  submitted ||
                  ['catalog_csv', 'catalog_tsv'].includes(inputFormat)
                }
              />
            )}
          />
          <Form.Check.Label className="font-weight-normal">
            WGS
          </Form.Check.Label>
        </Form.Check>
        <Form.Check inline id="radioWES">
          <Controller
            name="strategy"
            control={control}
            render={({ field }) => (
              <Form.Check.Input
                {...field}
                type="radio"
                value="WES"
                checked={strategy == 'WES'}
                disabled={
                  submitted ||
                  ['catalog_csv', 'catalog_tsv'].includes(inputFormat)
                }
              />
            )}
          />
          <Form.Check.Label className="font-weight-normal">
            WES
          </Form.Check.Label>
        </Form.Check>
      </Form.Group>
      <hr className="mb-3" />
      <Form.Group controlId="split" className="d-flex">
        <Form.Label className="mr-auto">
          Split Mutations According to Filter{' '}
          <OverlayTrigger
            trigger="click"
            placement="top"
            overlay={msPopover}
            rootClose
          >
            <Button
              variant="link"
              className="p-0 font-weight-bold"
              aria-label="mutation split info"
            >
              <FontAwesomeIcon
                icon={faInfoCircle}
                style={{ verticalAlign: 'baseline' }}
              />
            </Button>
          </OverlayTrigger>
        </Form.Label>
        <Form.Check inline id="split">
          <Controller
            name="mutationSplit"
            control={control}
            render={({ field }) => (
              <Form.Check.Input
                {...field}
                type="checkbox"
                disabled={
                  submitted ||
                  filter.length ||
                  bedFile.size ||
                  ['catalog_csv', 'catalog_tsv'].includes(inputFormat)
                }
              />
            )}
          />
        </Form.Check>
      </Form.Group>
      <Form.Group controlId="filter">
        <Form.Label>
          Select Filter{' '}
          <span className="text-muted font-italic font-weight-normal">
            (optional)
          </span>
        </Form.Label>
        <Controller
          name="filter"
          control={control}
          render={({ field }) => (
            <Form.Control
              {...field}
              type="text"
              size="sm"
              placeholder="Enter a filter"
              disabled={
                submitted ||
                mutationSplit ||
                ['catalog_csv', 'catalog_tsv'].includes(inputFormat)
              }
            />
          )}
        />
        <Form.Text className="text-muted">
          Use @ to separate multiple filters
        </Form.Text>
      </Form.Group>
      <hr className="mb-3" />
      <Form.Group controlId="dataFileUpload">
        <Form.Label>
          Filter Mutations using Bed File{' '}
          <span className="text-muted font-italic font-weight-normal">
            (optional)
          </span>
        </Form.Label>
        <Row className="m-0">
          <Col lg="6" className="p-0">
            <Button
              className="p-0 font-14"
              disabled={
                submitted ||
                mutationSplit ||
                ['catalog_csv', 'catalog_tsv'].includes(inputFormat)
              }
              variant="link"
              href={'assets/exampleInput/demo_input_bed.bed'}
              download
            >
              Download Example Bed
            </Button>
          </Col>
          <Col lg="6" className="p-0 d-flex">
            <Button
              className={`p-0 ml-auto font-14`}
              disabled={
                submitted ||
                isUploading ||
                isQueueing ||
                isLoading ||
                mutationSplit ||
                ['catalog_csv', 'catalog_tsv'].includes(inputFormat)
              }
              variant="link"
              type="button"
              onClick={() => loadBed()}
            >
              Load Example Bed
            </Button>
          </Col>
        </Row>
        <Row>
          <Col>
            <div className="d-flex">
              <Controller
                name="bedFile"
                control={control}
                render={({ field }) => (
                  <Form.File
                    disabled={
                      submitted ||
                      mutationSplit ||
                      ['catalog_csv', 'catalog_tsv'].includes(inputFormat) ||
                      isUploading ||
                      isQueueing ||
                      isLoading
                    }
                    id="uploadDataFile"
                    title={bedFilename || 'Upload Bed File...'}
                    value={''} // set dummy value for file input
                    label={
                      bedFile.size
                        ? bedFile.name
                        : bedFilename
                        ? bedFilename
                        : 'Upload Bed File...'
                    }
                    accept=".bed"
                    onChange={(e) => {
                      if (e.target.files.length) {
                        setValue('bedFile', e.target.files[0]);
                      }
                    }}
                    custom
                  />
                )}
              />
              {bedFile.size > 0 && (
                <Button
                  className="ml-1"
                  size="sm"
                  title="Remove"
                  variant="danger"
                  disabled={submitted || isUploading || isQueueing || isLoading}
                  onClick={removeBedFile}
                >
                  <FontAwesomeIcon icon={faFolderMinus} size="lg" />
                </Button>
              )}
            </div>
          </Col>
        </Row>
      </Form.Group>
      <hr className="mb-3" />
      <Form.Group controlId="collapse" className="d-flex">
        <Form.Label className="mr-auto">
          Add Collapsing Data{' '}
          <OverlayTrigger
            trigger="click"
            placement="top"
            overlay={csPopover}
            rootClose
          >
            <Button
              variant="link"
              className="p-0 font-weight-bold"
              aria-label="collapse data info"
            >
              <FontAwesomeIcon
                icon={faInfoCircle}
                style={{ verticalAlign: 'baseline' }}
              />
            </Button>
          </OverlayTrigger>
        </Form.Label>
        <Form.Check inline id="collapse">
          <Controller
            name="collapse"
            control={control}
            render={({ field }) => (
              <Form.Check.Input
                {...field}
                type="checkbox"
                disabled={
                  submitted ||
                  ['catalog_csv', 'catalog_tsv'].includes(inputFormat)
                }
              />
            )}
          />
        </Form.Check>
      </Form.Group>
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
            disabled={isUploading || isQueueing || isLoading}
            className="w-100"
            variant="secondary"
            onClick={() => handleReset()}
          >
            Reset
          </Button>
        </Col>
        <Col>
          <Button
            disabled={submitted || isUploading || isQueueing || isLoading}
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
