import React, { useState, useCallback } from 'react';
import { Form, Button, Row, Col } from 'react-bootstrap';
import { useDropzone } from 'react-dropzone';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faCloudUploadAlt, faMinus } from '@fortawesome/free-solid-svg-icons';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { useSelector } from 'react-redux';
import './visualize.scss';
import {
  dispatchVisualize,
  dispatchVisualizeResults,
  dispatchError,
  getInitialState,
} from '../../../services/store';
const { Group, Label, Control, Check, Text } = Form;

export default function UploadForm({ setOpenSidebar }) {
  const {
    inputFormat,
    selectedGenome,
    experimentalStrategy,
    mutationSplit,
    collapseSample,
    mutationFilter,
    queueMode,
    email,
    disableParameters,
    storeFile,
    submitted,
    exampleData,
    loading,
  } = useSelector((state) => state.visualize);

  const [inputFile, setInput] = useState(new File([], ''));

  const onDrop = useCallback((acceptedFiles) => {
    setInput(acceptedFiles[0]);
    dispatchVisualize({ storeFile: acceptedFiles[0].name });
  }, []);
  const { getRootProps, getInputProps } = useDropzone({
    onDrop,
    accept: '.csv, .tsv, .vcf, .gz, .zip, .tar, .tar.gz',
  });

  async function handleSubmit() {
    // disable parameters after submit
    dispatchVisualize({ disableParameters: true });
    dispatchVisualize({ submitted: true });

    const data = await uploadFile();
    if (data && data.projectID) {
      const args = {
        inputFormat: ['-f', inputFormat],
        inputFile: ['-i', data.filePath],
        projectID: ['-p', data.projectID],
        genomeAssemblyVersion: ['-g', selectedGenome],
        experimentalStrategy: ['-t', experimentalStrategy],
        collapseSample: ['-c', collapseSample],
        outputDir: ['-o', data.projectID],
      };
      // conditionally include params
      if (mutationFilter.length)
        args['mutationFilter'] = ['-F', mutationFilter];
      if (['vcf', 'csv', 'tsv'].includes(inputFormat))
        args['mutationSplit'] = ['-s', mutationSplit];

      if (queueMode) {
        dispatchVisualize({
          loading: {
            active: true,
            content: 'Sending to Queue...',
            showIndicator: true,
          },
        });
        try {
          const response = await fetch(`/visualize/queue`, {
            method: 'POST',
            headers: {
              Accept: 'application/json',
              'Content-Type': 'application/json',
            },
            body: JSON.stringify(data),
          });

          if (response.ok) {
            // placeholder alert with error modal
            dispatchError('Successfully submitted to queue.');
          } else {
            dispatchVisualizeResults({
              error: 'Please Reset Your Parameters and Try again.',
            });
            dispatchError('Failed to submit to queue. Please Try Again.');
          }
        } catch (err) {
          dispatchError(err);
        } finally {
          dispatchVisualize({ loading: { active: false } });
        }
      } else {
        dispatchVisualize({
          loading: {
            active: true,
            content: 'Calculating...',
            showIndicator: true,
          },
        });
        try {
          const response = await fetch(`/api/visualize`, {
            method: 'POST',
            headers: {
              Accept: 'application/json',
              'Content-Type': 'application/json',
            },
            body: JSON.stringify(args),
          });

          if (response.ok) {
            dispatchVisualizeResults({
              projectID: data.projectID,
              debug: await response.json(),
            });
            setOpenSidebar(false);
          } else if (response.status == 502) {
            dispatchVisualizeResults({
              error: 'Please Reset Your Parameters and Try again.',
            });
            dispatchError('Your submission has timed out. Please Try Again.');
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
        } finally {
          dispatchVisualize({
            loading: {
              active: false,
            },
          });
        }
      }
    }
  }

  function handleReset() {
    const initialState = getInitialState();
    // console.log(initialState.visualize);
    removeFile();
    dispatchVisualize(initialState.visualize);
    dispatchVisualizeResults(initialState.visualizeResults);
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
      data.append('file', inputFile);
      let response = await fetch(`/visualize/upload`, {
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
  }

  function selectFormat(format) {
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
    dispatchVisualize({ storeFile: filename });
  }

  return (
    <Form className="">
      <Group controlId="fileType">
        <Label>Choose File Type</Label>
        <Control
          as="select"
          value={inputFormat}
          onChange={(e) => selectFormat(e.target.value)}
          disabled={disableParameters}
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
        <Row className="m-0">
          <Col sm="6" className="p-0">
            <Button
              className="p-0"
              disabled={disableParameters}
              variant="link"
              href={exampleData}
              download
            >
              Download Sample
            </Button>
          </Col>
          <Col sm="6" className="p-0 d-flex">
            <Button
              className="p-0 ml-auto"
              disabled={disableParameters}
              variant="link"
              type="button"
              onClick={() => loadExample()}
            >
              Load Sample
            </Button>
          </Col>
        </Row>
        <section>
          <div {...getRootProps({ className: 'dropzone' })}>
            <input
              {...getInputProps()}
              disabled={inputFile.size || disableParameters}
            />
            {inputFile.size || submitted ? (
              <button
                id="removeFile"
                className="d-flex w-100 faButton"
                onClick={() => removeFile()}
                disabled={disableParameters}
              >
                <span id="uploadedFile">
                  {submitted ? storeFile : inputFile.name}
                </span>
                <span className="text-danger ml-auto">
                  <FontAwesomeIcon icon={faMinus} />
                </span>
              </button>
            ) : (
              <>
                <p>Drop files here or click to upload.</p>
                <FontAwesomeIcon icon={faCloudUploadAlt} size="4x" />
              </>
            )}
          </div>
        </section>
      </Group>
      <Group controlId="genomeAssembly">
        <Label>Genome Assembly Version</Label>
        <Control
          as="select"
          value={selectedGenome}
          onChange={(e) =>
            dispatchVisualize({ selectedGenome: e.target.value })
          }
          disabled={disableParameters}
          custom
        >
          <option value="GRCh37">GRCh37</option>
          <option value="GRCh38">GRCh38</option>
          <option value="mm10">mm10</option>
        </Control>
      </Group>
      <Group className="d-flex">
        <Label className="mr-auto">Experiment Strategy</Label>
        <Check inline id="radioWGS">
          <Check.Input
            type="radio"
            value="WGS"
            checked={experimentalStrategy == 'WGS'}
            onChange={(e) =>
              dispatchVisualize({ experimentalStrategy: e.target.value })
            }
            disabled={disableParameters}
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
            disabled={disableParameters}
          />
          <Check.Label className="font-weight-normal">WES</Check.Label>
        </Check>
      </Group>
      <Group className="d-flex">
        <Label className="mr-auto">Mutation Split</Label>
        <Check inline id="radioMutationSplitFalse">
          <Check.Input
            disabled={!['vcf', 'csv', 'tsv'].includes(inputFormat)}
            type="radio"
            value="False"
            checked={mutationSplit == 'False'}
            onChange={(e) =>
              dispatchVisualize({ mutationSplit: e.target.value })
            }
            disabled={disableParameters}
          />
          <Check.Label className="font-weight-normal">False</Check.Label>
        </Check>
        <Check inline id="radioMutationSplitTrue">
          <Check.Input
            disabled={!['vcf', 'csv', 'tsv'].includes(inputFormat)}
            type="radio"
            value="True"
            checked={mutationSplit == 'True'}
            onChange={(e) =>
              dispatchVisualize({ mutationSplit: e.target.value })
            }
            disabled={disableParameters}
          />
          <Check.Label className="font-weight-normal">True</Check.Label>
        </Check>
      </Group>

      <Group className="d-flex">
        <Label className="mr-auto">Collapse Sample</Label>
        <Check inline id="radioFalse">
          <Check.Input
            disabled={disableParameters}
            type="radio"
            value="False"
            checked={collapseSample == 'False'}
            onChange={(e) =>
              dispatchVisualize({ collapseSample: e.target.value })
            }
          />
          <Check.Label className="font-weight-normal">False</Check.Label>
        </Check>
        <Check inline id="radioTrue">
          <Check.Input
            disabled={disableParameters}
            type="radio"
            value="True"
            checked={collapseSample == 'True'}
            onChange={(e) =>
              dispatchVisualize({ collapseSample: e.target.value })
            }
          />
          <Check.Label className="font-weight-normal">True</Check.Label>
        </Check>
      </Group>
      <Group controlId="filter">
        <Label>Mutation Filter</Label>
        <Control
          type="text"
          size="sm"
          placeholder="Enter a filter"
          value={mutationFilter}
          onChange={(e) =>
            dispatchVisualize({ mutationFilter: e.target.value })
          }
          disabled={disableParameters}
        ></Control>
        <Text className="text-muted">Use @ to separate multiple filters</Text>
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
            disabled={!inputFile.size || disableParameters || loading.active}
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
