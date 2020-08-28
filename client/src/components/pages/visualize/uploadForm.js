import React, { useState, useCallback } from 'react';
import {
  Form,
  Button,
  Row,
  Col,
  Popover,
  OverlayTrigger,
} from 'react-bootstrap';
import { useDropzone } from 'react-dropzone';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import {
  faCloudUploadAlt,
  faTimes,
  faInfoCircle,
} from '@fortawesome/free-solid-svg-icons';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { useSelector } from 'react-redux';
import './visualize.scss';
import {
  dispatchVisualize,
  dispatchVisualizeResults,
  dispatchError,
  getInitialState,
  dispatchMutationalProfiles,
  dispatchCosineSimilarity,
  dispatchProfileComparison,
  dispatchPCA,
} from '../../../services/store';
const { Group, Label, Control, Check, Text } = Form;
const { Title, Content } = Popover;

export default function UploadForm() {
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
  } = useSelector((state) => state.visualize);
  const rootURL = window.location.pathname;
  const [inputFile, setInput] = useState(new File([], ''));
  const [bedFile, setBed] = useState(new File([], ''));
  const onDropMain = useCallback((acceptedFiles) => {
    setInput(acceptedFiles[0]);
    dispatchVisualize({ storeFilename: acceptedFiles[0].name });
  }, []);
  const {
    getRootProps: mainRootProps,
    getInputProps: mainInputProps,
  } = useDropzone({
    onDrop: onDropMain,
    accept: '.csv, .tsv, .vcf, .gz, .zip, .tar, .tar.gz',
  });
  const onDropBed = useCallback((acceptedFiles) => {
    setBed(acceptedFiles[0]);
    dispatchVisualize({ bedFilename: acceptedFiles[0].name });
  }, []);
  const {
    getRootProps: bedRootProps,
    getInputProps: bedInputProps,
  } = useDropzone({
    onDrop: onDropBed,
    accept: '.bed',
  });

  async function handleSubmit() {
    // disable parameters after submit
    dispatchVisualize({ submitted: true });

    const { projectID, filePath, bedPath } = await uploadFile();
    if (projectID) {
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
        // wip
        // dispatchVisualize({
        //   loading: {
        //     active: true,
        //     content: 'Sending to Queue...',
        //     showIndicator: true,
        //   },
        // });
        // try {
        //   const response = await fetch(`${rootURL}visualize/queue`, {
        //     method: 'POST',
        //     headers: {
        //       Accept: 'application/json',
        //       'Content-Type': 'application/json',
        //     },
        //     body: JSON.stringify(data),
        //   });
        //   if (response.ok) {
        //     // placeholder alert with error modal
        //     dispatchError('Successfully submitted to queue.');
        //   } else {
        //     dispatchVisualizeResults({
        //       error: 'Please Reset Your Parameters and Try again.',
        //     });
        //     dispatchError('Failed to submit to queue. Please Try Again.');
        //   }
        // } catch (err) {
        //   dispatchError(err);
        // }
      } else {
        dispatchVisualize({
          loading: {
            active: true,
            content: 'Calculating...',
            showIndicator: true,
          },
        });
        try {
          const response = await fetch(`${rootURL}api/visualize`, {
            method: 'POST',
            headers: {
              Accept: 'application/json',
              'Content-Type': 'application/json',
            },
            body: JSON.stringify(args),
          });

          if (response.ok) {
            const results = await response.json();
            dispatchVisualizeResults({
              projectID: projectID,
              svgList: results.svgList,
              statistics: results.statistics,
              matrixList: results.matrixList,
              downloads: results.downloads,
            });
            dispatchMutationalProfiles({
              debug: { stdout: results.stdout, stderr: results.stderr },
            });
          } else if (response.status == 504) {
            dispatchVisualizeResults({
              error: 'Please Reset Your Parameters and Try again.',
            });
            dispatchError(
              'Your submission has timed out. Please try again by submitting this job to a queue instead.'
            );
          } else {
            dispatchVisualizeResults({
              error: 'Please Reset Your Parameters and Try again.',
            });
            const { msg, stdout, stderr } = await response.json();
            const message = `<div>
            <p>${msg}</p>
            <p><b>Python:</b></p>
            <pre>${stdout}</pre>
            <pre>${stderr}</pre>
          </div>`;
            dispatchError(message);
          }
        } catch (err) {
          dispatchError(err);
          dispatchVisualizeResults({
            error: 'Please Reset Your Parameters and Try again.',
          });
        }
        dispatchVisualize({ loading: { active: false } });
      }
    }
  }

  function handleReset() {
    const initialState = getInitialState();
    removeFile();
    removeBedFile();
    dispatchVisualize(initialState.visualize);
    dispatchVisualizeResults(initialState.visualizeResults);
    dispatchMutationalProfiles(initialState.mutationalProfiles);
    dispatchCosineSimilarity(initialState.cosineSimilarity);
    dispatchProfileComparison(initialState.profileComparison);
    dispatchPCA(initialState.pca);
  }

  //   Uploads inputFile and returns a projectID
  async function uploadFile() {
    dispatchVisualize({
      loading: {
        active: true,
        content: 'Uploading file...',
        showIndicator: true,
      },
    });
    try {
      const data = new FormData();
      data.append('inputFile', inputFile);
      if (bedFile.size) data.append('bedFile', bedFile);
      let response = await fetch(`${rootURL}visualize/upload`, {
        method: 'POST',
        body: data,
      });

      if (!response.ok) {
        const { msg, error } = await response.json();
        const message = `<div>
          <p>${msg}</p>
         ${error ? `<p>${error}</p>` : ''} 
        </div>`;
        dispatchError(message);
      } else {
        return await response.json();
      }
    } catch (err) {
      dispatchError(err);
    } finally {
      dispatchVisualize({
        loading: {
          active: false,
        },
      });
    }
  }

  function removeFile() {
    setInput(new File([], ''));
    dispatchVisualize({ storeFilename: '' });
  }

  function removeBedFile() {
    setBed(new File([], ''));
    dispatchVisualize({ bedFilename: '' });
  }

  function selectFormat(format) {
    handleReset();
    let path = '';
    if (format == 'vcf') path = 'assets/exampleInput/demo_input_multi.vcf.gz';
    if (format == 'csv') path = 'assets/exampleInput/demo_input_multi.csv';
    if (format == 'tsv') path = 'assets/exampleInput/demo_input_multi.tsv';
    if (format == 'catalog_tsv')
      path = 'assets/exampleInput/demo_input_catalog.tsv';
    if (format == 'catalog_csv')
      path = 'assets/exampleInput/demo_input_catalog.csv';

    dispatchVisualize({ inputFormat: format, exampleData: path });
  }

  async function loadExample() {
    const filename = exampleData.split('/').slice(-1)[0];
    setInput(new File([await (await fetch(exampleData)).blob()], filename));
    dispatchVisualize({ storeFilename: filename });
  }

  async function loadBed() {
    const filename = bedData.split('/').slice(-1)[0];
    setBed(new File([await (await fetch(bedData)).blob()], filename));
    dispatchVisualize({ bedFilename: filename });
  }

  const msPopover = (
    <Popover id="popover-basic">
      <Title as="h3">Mutation Split</Title>
      <Content>
        <p>
          For each sample, split mutations in to different groups according the
          “Filter” column in VCF/CSV/TSV file. Splitting operation use the “;”
          as separator.
        </p>
      </Content>
    </Popover>
  );

  const csPopover = (
    <Popover id="popover-basic">
      <Title as="h3">Mutation Split</Title>
      <Content>
        <p>
          A new sample called “All_Sample” will add to the result, which
          combines the mutations from all samples.
        </p>
      </Content>
    </Popover>
  );

  return (
    <Form className="">
      <Group controlId="fileType">
        <Label>Choose File Format</Label>
        <Control
          as="select"
          value={inputFormat}
          onChange={(e) => selectFormat(e.target.value)}
          disabled={submitted}
          custom
        >
          <option value="vcf">VCF</option>
          <option value="csv">CSV</option>
          <option value="tsv">TSV</option>
          <option value="catalog_tsv">CATALOG TSV</option>
          <option value="catalog_csv">CATALOG CSV</option>
        </Control>
      </Group>
      <Group controlId="fileUpload">
        <Label>
          Upload File <span style={{ color: 'red' }}>*</span>
        </Label>
        <Row className="m-0">
          <Col sm="6" className="p-0">
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
          <Col sm="6" className="p-0 d-flex">
            <Button
              className="p-0 ml-auto font-14"
              disabled={submitted}
              variant="link"
              type="button"
              onClick={() => loadExample()}
            >
              Load Example Data
            </Button>
          </Col>
        </Row>
        <section>
          <div
            {...mainRootProps({ className: 'dropzone' })}
            disabled={submitted}
          >
            <input
              {...mainInputProps()}
              disabled={inputFile.size || submitted}
            />
            {inputFile.size ? (
              <button
                id="removeFile"
                className="d-flex w-100 faButton"
                onClick={() => removeFile()}
                disabled={submitted}
              >
                <span id="uploadedFile">
                  {submitted ? storeFilename : inputFile.name}
                </span>
                <span className="text-danger ml-auto">
                  <FontAwesomeIcon icon={faTimes} />
                </span>
              </button>
            ) : (
              <>
                <span>Drop files here or click to upload</span>
                {/* <FontAwesomeIcon icon={faCloudUploadAlt} size="4x" /> */}
              </>
            )}
          </div>
        </section>
      </Group>
      <Group controlId="genomeAssembly">
        <Label>Choose Reference Genome Build</Label>
        <Control
          as="select"
          value={selectedGenome}
          onChange={(e) =>
            dispatchVisualize({ selectedGenome: e.target.value })
          }
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
              dispatchVisualize({ experimentalStrategy: e.target.value })
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
              dispatchVisualize({ experimentalStrategy: e.target.value })
            }
            disabled={
              submitted || ['catalog_csv', 'catalog_tsv'].includes(inputFormat)
            }
          />
          <Check.Label className="font-weight-normal">WES</Check.Label>
        </Check>
      </Group>
      <hr />
      <Group className="d-flex">
        <Label className="mr-auto">
          Split Mutations According to Filter{' '}
          <OverlayTrigger
            trigger="click"
            placement="top"
            overlay={msPopover}
            rootClose
          >
            <Button variant="link" className="p-0 font-weight-bold">
              <FontAwesomeIcon
                icon={faInfoCircle}
                style={{ verticalAlign: 'baseline' }}
              />
            </Button>
          </OverlayTrigger>
        </Label>
        <Check inline id="radioMutationSplitFalse">
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
              dispatchVisualize({
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
          onChange={(e) =>
            dispatchVisualize({ mutationFilter: e.target.value })
          }
          disabled={
            submitted ||
            mutationSplit == 'True' ||
            ['catalog_csv', 'catalog_tsv'].includes(inputFormat)
          }
        ></Control>
        <Text className="text-muted">Use @ to separate multiple filters</Text>
      </Group>
      <hr />
      <Group controlId="bedUpload">
        <Label>
          Filter Mutations using Bed File{' '}
          <span className="text-muted font-italic font-weight-normal">
            (optional)
          </span>
        </Label>
        <Row className="m-0">
          <Col sm="6" className="p-0">
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
              Download Example Bed Data
            </Button>
          </Col>
          <Col sm="6" className="p-0 d-flex">
            <Button
              className="p-0 ml-auto font-14"
              disabled={
                submitted ||
                mutationSplit == 'True' ||
                ['catalog_csv', 'catalog_tsv'].includes(inputFormat)
              }
              variant="link"
              type="button"
              onClick={() => loadBed()}
            >
              Load Example Bed Data
            </Button>
          </Col>
        </Row>
        <section>
          <div
            disabled={
              submitted ||
              mutationSplit == 'True' ||
              ['catalog_csv', 'catalog_tsv'].includes(inputFormat)
            }
            {...bedRootProps({ className: 'dropzone' })}
          >
            <input
              {...bedInputProps()}
              disabled={
                bedFile.size ||
                submitted ||
                ['catalog_csv', 'catalog_tsv'].includes(inputFormat)
              }
            />
            {bedFile.size ? (
              bedFilename.length > 0 && (
                <button
                  id="removeFile"
                  className="d-flex w-100 faButton"
                  onClick={() => removeBedFile()}
                  disabled={submitted}
                >
                  <span id="uploadedFile">
                    {submitted ? bedFilename : bedFile.name}
                  </span>
                  <span className="text-danger ml-auto">
                    <FontAwesomeIcon icon={faTimes} />
                  </span>
                </button>
              )
            ) : (
              <>
                <span>Drop files here or click to upload</span>
                {/* <FontAwesomeIcon icon={faCloudUploadAlt} size="4x" /> */}
              </>
            )}
          </div>
        </section>
      </Group>
      <hr />
      <Group className="d-flex">
        <Label className="mr-auto">
          Add Collapsing Data{' '}
          <OverlayTrigger
            trigger="click"
            placement="top"
            overlay={csPopover}
            rootClose
          >
            <Button variant="link" className="p-0 font-weight-bold">
              <FontAwesomeIcon
                icon={faInfoCircle}
                style={{ verticalAlign: 'baseline' }}
              />
            </Button>
          </OverlayTrigger>
        </Label>
        <Check inline id="radioFalse">
          <Check.Input
            disabled={
              submitted || ['catalog_csv', 'catalog_tsv'].includes(inputFormat)
            }
            type="checkbox"
            value={collapseSample}
            checked={collapseSample == 'True'}
            onChange={(e) =>
              dispatchVisualize({
                collapseSample: e.target.value == 'True' ? 'False' : 'True',
              })
            }
          />
        </Check>
      </Group>
      <hr />
      <Group controlId="email">
        <LoadingOverlay active={true} content={'Work in progress...'} />
        <div>
          <Check
            inline
            id="toggleQueue"
            type="checkbox"
            label="Submit this job to a Queue"
            checked={queueMode == true}
            onChange={(_) => {
              dispatchVisualize({ queueMode: !queueMode });
            }}
          />
        </div>
        <div>
          <Control
            placeholder="Enter Email"
            size="sm"
            value={email}
            onChange={(e) => dispatchVisualize({ email: e.target.value })}
            disabled={!queueMode}
          ></Control>
          <Text className="text-muted">
            <i>
              Note: if sending to queue, when computation is completed, a
              notification will be sent to the e-mail entered above.
            </i>
          </Text>
        </div>
      </Group>
      <Row>
        <Col sm="6">
          <Button
            disabled={loading.active}
            className="w-100"
            variant="secondary"
            onClick={(e) => handleReset()}
          >
            Reset
          </Button>
        </Col>
        <Col sm="6">
          <Button
            disabled={!inputFile.size || submitted || loading.active}
            className="w-100"
            variant="primary"
            type="button"
            onClick={(e) => handleSubmit(e)}
          >
            Submit
          </Button>
        </Col>
      </Row>
    </Form>
  );
}
