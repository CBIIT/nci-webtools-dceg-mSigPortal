import React, { useState } from 'react';
import {
  Form,
  Button,
  Row,
  Col,
  Popover,
  OverlayTrigger,
} from 'react-bootstrap';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faInfoCircle, faFolderMinus } from '@fortawesome/free-solid-svg-icons';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import './visualization.scss';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../services/store/visualization';
import { actions as modalActions } from '../../../services/store/modal';

const actions = { ...visualizationActions, ...modalActions };
const { Group, Label, Control, Check, Text } = Form;
const { Title, Content } = Popover;

export default function UserForm() {
  const dispatch = useDispatch();
  const visualization = useSelector((state) => state.visualization);
  const mergeState = (state) =>
    dispatch(actions.mergeVisualization({ state: state }));
  const mergeMutationalProfiles = (state) =>
    dispatch(actions.mergeVisualization({ mutationalProfiles: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));
  const mergeSuccess = (msg) =>
    dispatch(actions.mergeModal({ success: { visible: true, message: msg } }));
  const resetVisualization = (_) => dispatch(actions.resetVisualization());

  const {
    inputFormat,
    selectedGenome,
    experimentalStrategy,
    mutationSplit,
    collapseSample,
    mutationFilter,
    queueMode,
    email,
    storeFilename,
    bedFilename,
    bedData,
    submitted,
    exampleData,
    loading,
    pDataOptions,
    studyOptions,
    study,
    cancerTypeOptions,
    cancerType,
    pubExperimentOptions,
    pubExperimentalStrategy,
  } = visualization.main;

  const [inputFile, setInput] = useState(new File([], ''));
  const [bedFile, setBed] = useState(new File([], ''));
  const [validFile, setValidFile] = useState(false);
  const [validEmail, setValidEmail] = useState(false);
  const [checkValid, setCheckValid] = useState(false);

  async function handleSubmit() {
    const { projectID, filePath, bedPath } = await uploadFile();

    const args = {
      inputFormat: ['-f', inputFormat],
      inputFile: ['-i', filePath],
      projectID: ['-p', projectID],
      genomeAssemblyVersion: ['-g', selectedGenome],
      experimentalStrategy: ['-t', experimentalStrategy],
      outputDir: ['-o', projectID],
    };

    // conditionally include mutation split and mutation filter params
    if (['vcf', 'csv', 'tsv'].includes(inputFormat)) {
      args['collapseSample'] = ['-c', collapseSample];

      if (bedFile.size) args['bedFile'] = ['-b', bedPath];

      if (mutationFilter.length) {
        args['mutationFilter'] = ['-F', mutationFilter];
      } else {
        args['mutationSplit'] = ['-s', mutationSplit];
      }
    }

    if (queueMode) {
      mergeState({
        loading: {
          active: true,
          // content: 'Sending to Queue...',
          // showIndicator: true,
        },
        submitted: true,
      });

      try {
        const response = await fetch(`web/queue`, {
          method: 'POST',
          headers: {
            Accept: 'application/json',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            args: args,
            state: {
              ...visualization,
            },
          }),
        });

        mergeState({ loading: { active: false } });

        if (response.ok) {
          // placeholder alert with error modal
          mergeSuccess(
            `Your job was successfully submitted to the queue. You will recieve an email at ${email} with your results.`
          );
        } else {
          mergeState({
            error: 'Please Reset Your Parameters and Try again.',
            submitted: false,
          });
          mergeError('Failed to submit to queue. Please Try Again.');
        }
      } catch (err) {
        mergeError(err.message);
        mergeState({ loading: { active: false } });
      }
    } else {
      mergeState({
        loading: {
          active: true,
          // content: 'Calculating...',
          // showIndicator: true,
        },
      });
      try {
        const response = await fetch(`web/profilerExtraction`, {
          method: 'POST',
          headers: {
            Accept: 'application/json',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify(args),
        });

        if (response.ok) {
          const results = await response.json();

          mergeState({
            projectID: projectID,
            svgList: results.svgList,
            statistics: results.statistics,
            matrixList: results.matrixList,
            downloads: results.downloads,
          });
          mergeMutationalProfiles({
            debug: { stdout: results.stdout, stderr: results.stderr },
          });
        } else if (response.status == 504) {
          mergeState({
            error: 'Please Reset Your Parameters and Try again.',
          });
          mergeError(
            'Your submission has timed out. Please try again by submitting this job to a queue instead.'
          );
        } else {
          mergeState({
            error: 'Please Reset Your Parameters and Try again.',
          });
          const { stdout, stderr } = await response.json();

          const message = `<div>
            <pre>${stdout}</pre>
            <pre>${stderr}</pre>
          </div>`;
          mergeError(message);
        }
      } catch (err) {
        mergeError(err.message);
        mergeState({
          error: 'Please Reset Your Parameters and Try again.',
        });
      }
      mergeState({
        loading: { active: false },
        submitted: true,
        displayTab: 'profilerSummary',
        openSidebar: false,
      });
    }
  }

  // reset form while preserving selected data source and public parameters
  function handleReset() {
    const params = {
      source: 'user',
      pDataOptions: pDataOptions,
      studyOptions: studyOptions,
      study: study,
      cancerTypeOptions: cancerTypeOptions,
      cancerType: cancerType,
      pubExperimentOptions: pubExperimentOptions,
      pubExperimentalStrategy: pubExperimentalStrategy,
    };
    // clear id from url
    window.location.hash = '#/visualization';
    resetVisualization();
    setCheckValid(false);
    removeFile();
    removeBedFile();
    mergeState(params);
  }

  //   Uploads inputFile and returns a projectID
  async function uploadFile() {
    mergeState({
      loading: {
        active: true,
        // content: 'Uploading file...',
        // showIndicator: true,
      },
    });
    try {
      const data = new FormData();
      data.append('inputFile', inputFile);
      if (bedFile.size) data.append('bedFile', bedFile);
      let response = await fetch(`web/upload`, {
        method: 'POST',
        body: data,
      });

      if (!response.ok) {
        const { msg, error } = await response.json();
        const message = `<div>
          <p>${msg}</p>
         ${error ? `<p>${error}</p>` : ''} 
        </div>`;
        mergeError(message);
      } else {
        return await response.json();
      }
    } catch (err) {
      mergeError(err.message);
    } finally {
      mergeState({
        loading: {
          active: false,
        },
      });
    }
  }

  function removeFile() {
    setInput(new File([], ''));
    mergeState({ storeFilename: '' });
  }

  function removeBedFile() {
    setBed(new File([], ''));
    mergeState({ bedFilename: '' });
  }

  function selectFormat(format) {
    removeFile();
    removeBedFile();
    let path = '';
    if (format == 'vcf') path = 'assets/exampleInput/demo_input_multi.vcf.gz';
    if (format == 'maf') path = 'assets/exampleInput/demo_input_multi_MAF.txt';
    if (format == 'csv') path = 'assets/exampleInput/demo_input_multi.csv';
    if (format == 'tsv') path = 'assets/exampleInput/demo_input_multi.tsv';
    if (format == 'catalog_tsv')
      path = 'assets/exampleInput/demo_input_catalog.tsv';
    if (format == 'catalog_csv')
      path = 'assets/exampleInput/demo_input_catalog.csv';

    mergeState({ inputFormat: format, exampleData: path });
  }

  async function loadExample() {
    const filename = exampleData.split('/').slice(-1)[0];
    if (storeFilename != filename) {
      setInput(new File([await (await fetch(exampleData)).blob()], filename));
      mergeState({ storeFilename: filename });
    }
  }

  async function loadBed() {
    const filename = bedData.split('/').slice(-1)[0];
    if (bedFilename != filename) {
      setBed(new File([await (await fetch(bedData)).blob()], filename));
      mergeState({ bedFilename: filename });
    }
  }

  function validateForm() {
    setCheckValid(true);
    const re = new RegExp(
      /^[a-zA-Z0-9.!#$%&'*+/=?^_`{|}~-]+@[a-zA-Z0-9-]+(?:\.[a-zA-Z0-9-]+)*$/
    );
    setValidFile(inputFile.size > 0);
    if (queueMode) {
      setValidEmail(re.test(email));
      return inputFile.size > 0 && re.test(email);
    } else {
      return inputFile.size > 0;
    }
  }

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
    <Form>
      <Group controlId="fileType">
        <Label>Select File Format</Label>
        <Control
          as="select"
          value={inputFormat}
          onChange={(e) => selectFormat(e.target.value)}
          disabled={submitted}
          custom
        >
          <option value="vcf">VCF</option>
          <option value="maf">MAF</option>
          <option value="csv">CSV</option>
          <option value="tsv">TSV</option>
          <option value="catalog_tsv">CATALOG TSV</option>
          <option value="catalog_csv">CATALOG CSV</option>
        </Control>
      </Group>
      <Group>
        <Label>
          Upload File <span style={{ color: 'crimson' }}>*</span>
        </Label>
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
              disabled={submitted || loading.active}
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
              <Form.File
                disabled={submitted || loading.active}
                id="fileUpload"
                title={storeFilename || 'Upload Data File...'}
                value={''}
                label={
                  inputFile.size
                    ? inputFile.name
                    : storeFilename
                    ? storeFilename
                    : 'Upload Data File...'
                }
                accept=".csv, .tsv, .vcf, .gz, .tar, .tar.gz, .txt"
                isInvalid={checkValid ? !validFile : false}
                feedback="Please upload a data file"
                onChange={(e) => {
                  if (e.target.files.length) {
                    setInput(e.target.files[0]);
                    mergeState({
                      storeFilename: e.target.files[0].name,
                    });
                  }
                }}
                custom
              />
              {inputFile.size > 0 && (
                <Button
                  className="ml-1"
                  size="sm"
                  title="Remove"
                  variant="danger"
                  disabled={submitted || loading.active}
                  onClick={removeFile}
                >
                  <FontAwesomeIcon icon={faFolderMinus} size="lg" />
                </Button>
              )}
            </div>
          </Col>
        </Row>
      </Group>
      <Group controlId="genomeAssembly">
        <Label>Select Reference Genome Build</Label>
        <Control
          as="select"
          value={selectedGenome}
          onChange={(e) => mergeState({ selectedGenome: e.target.value })}
          disabled={
            submitted || ['catalog_csv', 'catalog_tsv'].includes(inputFormat)
          }
          custom
        >
          <option value="GRCh37">GRCh37</option>
          <option value="GRCh38">GRCh38</option>
          <option value="mm10">mm10</option>
        </Control>
      </Group>
      <Group className="d-flex">
        <Label className="mr-auto">Experimental Strategy</Label>
        <Check inline id="radioWGS">
          <Check.Input
            type="radio"
            value="WGS"
            checked={experimentalStrategy == 'WGS'}
            onChange={(e) =>
              mergeState({ experimentalStrategy: e.target.value })
            }
            disabled={
              submitted || ['catalog_csv', 'catalog_tsv'].includes(inputFormat)
            }
          />
          <Check.Label className="font-weight-normal">WGS</Check.Label>
        </Check>
        <Check inline id="radioWES">
          <Check.Input
            type="radio"
            value="WES"
            checked={experimentalStrategy == 'WES'}
            onChange={(e) =>
              mergeState({ experimentalStrategy: e.target.value })
            }
            disabled={
              submitted || ['catalog_csv', 'catalog_tsv'].includes(inputFormat)
            }
          />
          <Check.Label className="font-weight-normal">WES</Check.Label>
        </Check>
      </Group>
      <hr className="mb-3" />
      <Group controlId="split" className="d-flex">
        <Label className="mr-auto">
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
        </Label>
        <Check inline id="split">
          <Check.Input
            disabled={
              submitted ||
              mutationFilter.length ||
              bedFile.size ||
              ['catalog_csv', 'catalog_tsv'].includes(inputFormat)
            }
            type="checkbox"
            value={mutationSplit}
            checked={mutationSplit == 'True'}
            onChange={(e) =>
              mergeState({
                mutationSplit: e.target.value == 'True' ? 'False' : 'True',
              })
            }
          />
        </Check>
      </Group>
      <Group controlId="filter">
        <Label>
          Select Filter{' '}
          <span className="text-muted font-italic font-weight-normal">
            (optional)
          </span>
        </Label>
        <Control
          type="text"
          size="sm"
          placeholder="Enter a filter"
          value={mutationFilter}
          onChange={(e) => mergeState({ mutationFilter: e.target.value })}
          disabled={
            submitted ||
            mutationSplit == 'True' ||
            ['catalog_csv', 'catalog_tsv'].includes(inputFormat)
          }
        ></Control>
        <Text className="text-muted">Use @ to separate multiple filters</Text>
      </Group>
      <hr className="mb-3" />
      <Group controlId="dataFileUpload">
        <Label>
          Filter Mutations using Bed File{' '}
          <span className="text-muted font-italic font-weight-normal">
            (optional)
          </span>
        </Label>
        <Row className="m-0">
          <Col lg="6" className="p-0">
            <Button
              className="p-0 font-14"
              disabled={
                submitted ||
                mutationSplit == 'True' ||
                ['catalog_csv', 'catalog_tsv'].includes(inputFormat)
              }
              variant="link"
              href={bedData}
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
                loading.active ||
                mutationSplit == 'True' ||
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
              <Form.File
                disabled={
                  submitted ||
                  mutationSplit == 'True' ||
                  ['catalog_csv', 'catalog_tsv'].includes(inputFormat) ||
                  loading.active
                }
                id="uploadDataFile"
                title={bedFilename || 'Upload Bed File...'}
                value={''}
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
                    setBed(e.target.files[0]);
                    mergeState({
                      bedFilename: e.target.files[0].name,
                    });
                  }
                }}
                custom
              />
              {bedFile.size > 0 && (
                <Button
                  className="ml-1"
                  size="sm"
                  title="Remove"
                  variant="danger"
                  disabled={submitted || loading.active}
                  onClick={removeBedFile}
                >
                  <FontAwesomeIcon icon={faFolderMinus} size="lg" />
                </Button>
              )}
            </div>
          </Col>
        </Row>
      </Group>
      <hr className="mb-3" />
      <Group controlId="collapse" className="d-flex">
        <Label className="mr-auto">
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
        </Label>
        <Check inline id="collapse">
          <Check.Input
            disabled={
              submitted || ['catalog_csv', 'catalog_tsv'].includes(inputFormat)
            }
            type="checkbox"
            value={collapseSample}
            checked={collapseSample == 'True'}
            onChange={(e) =>
              mergeState({
                collapseSample: e.target.value == 'True' ? 'False' : 'True',
              })
            }
          />
        </Check>
      </Group>
      <hr className="mb-3" />
      <div>
        <LoadingOverlay active={false} content={'Work in progress...'} />
        <Group controlId="toggleQueue" className="d-flex">
          <Label className="mr-auto">Submit this job to a Queue</Label>{' '}
          <Check inline>
            <Check.Input
              type="checkbox"
              disabled={submitted}
              checked={queueMode}
              onChange={(_) => {
                mergeState({ queueMode: !queueMode });
              }}
            />
          </Check>
        </Group>
        <Group controlId="email">
          <Control
            aria-label="Enter Email"
            placeholder="Enter Email"
            size="sm"
            value={email}
            type="email"
            onChange={(e) => mergeState({ email: e.target.value })}
            disabled={!queueMode || submitted}
            isInvalid={queueMode && checkValid ? !validEmail : false}
          />
          <Control.Feedback type="invalid">
            Please provide a valid email
          </Control.Feedback>
          <Text className="text-muted">
            <i>
              Note: If sending to queue, when computation is completed, a
              notification will be sent to the e-mail entered above.
            </i>
          </Text>
        </Group>
      </div>
      <Row>
        <Col>
          <Button
            disabled={loading.active}
            className="w-100"
            variant="secondary"
            onClick={() => handleReset()}
          >
            Reset
          </Button>
        </Col>
        <Col>
          <Button
            disabled={submitted || loading.active}
            className="w-100"
            variant="primary"
            type="button"
            onClick={() => {
              if (validateForm()) handleSubmit();
            }}
          >
            Submit
          </Button>
        </Col>
      </Row>
    </Form>
  );
}
