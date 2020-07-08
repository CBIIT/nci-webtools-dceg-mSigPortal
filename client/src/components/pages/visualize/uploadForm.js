import React, { useState, useCallback } from 'react';
import { Form, Button } from 'react-bootstrap';
import { useDropzone } from 'react-dropzone';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faCloudUploadAlt, faMinus } from '@fortawesome/free-solid-svg-icons';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { useSelector } from 'react-redux';
import './visualize.scss';
import {
  store,
  updateVisualize,
  updateVisualizeResults,
  updateError,
  getInitialState,
} from '../../../services/store';
const { Group, Label, Control, Check, Text } = Form;

const root =
  process.env.NODE_ENV === 'development'
    ? 'http://localhost:8330/'
    : window.location.pathname;

export default function UploadForm({ setOpenSidebar }) {
  const {
    inputFormat,
    selectedGenome,
    experimentalStrategy,
    mutationSplit,
    isMultiple,
    collapseSample,
    mutationFilter,
    queueMode,
    email,
    disableParameters,
    storeFile,
    loading,
  } = useSelector((state) => state.visualize);

  const [inputFile, setInput] = useState(new File([], ''));

  const onDrop = useCallback((acceptedFiles) => {
    setInput(acceptedFiles[0]);
    store.dispatch(updateVisualize({ storeFile: acceptedFiles[0].name }));
  }, []);
  const { getRootProps, getInputProps } = useDropzone({
    onDrop,
    accept: '.csv, .tsv, .vcf, .gz, .zip, .tar, .tar.gz',
  });

  async function handleSubmit() {
    // disable parameters after submit
    store.dispatch(updateVisualize({ disableParameters: true }));

    const data = await uploadFile();
    if (data && data.projectID) {
      const args = {
        inputFormat: ['-f', inputFormat],
        inputFile: ['-i', data.filePath],
        projectID: ['-p', data.projectID],
        genomeAssemblyVersion: ['-g', selectedGenome],
        experimentalStrategy: ['-t', experimentalStrategy],
        outputDir: ['-o', data.projectID],
      };
      // conditionally include params
      if (mutationFilter.length)
        args['mutationFilter'] = ['-F', mutationFilter];
      if (['vcf', 'csv', 'tsv'].includes(inputFormat))
        args['mutationSplit'] = ['-s', mutationSplit];
      if (isMultiple) args['collapseSample'] = ['-c', collapseSample];

      if (queueMode) {
        store.dispatch(
          updateVisualize({
            loading: {
              active: true,
              content: 'Sending to Queue...',
              showIndicator: true,
            },
          })
        );
      } else {
        store.dispatch(
          updateVisualize({
            loading: {
              active: true,
              content: 'Calculating...',
              showIndicator: true,
            },
          })
        );

        const response = await fetch(`${root}api/visualize`, {
          method: 'POST',
          headers: {
            Accept: 'application/json',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify(args),
        });

        if (response.ok) {
          store.dispatch(
            updateVisualizeResults({
              projectID: data.projectID,
            })
          );
          setOpenSidebar(false);
        } else if (response.status == 502) {
          store.dispatch(
            updateVisualizeResults({
              error: 'Please Reset Your Parameters and Try again.',
            })
          );
          store.dispatch(
            updateError({
              visible: true,
              message: 'Your submission has timed out. Please Try Again.',
            })
          );
        } else {
          store.dispatch(
            updateVisualizeResults({
              error: 'Please Reset Your Parameters and Try again.',
            })
          );
          const { msg, stdout, stderr } = await response.json();
          const message = `<div>
            <p>${msg}</p>
            <p><b>Python:</b></p>
            <pre>${stdout}</pre>
            <pre>${stderr}</pre>
          </div>`;
          store.dispatch(updateError({ visible: true, message: message }));
        }
        store.dispatch(
          updateVisualize({
            loading: {
              active: false,
              showIndicator: false,
            },
          })
        );
      }
    }
  }

  function handleReset() {
    const initialState = getInitialState();
    // console.log(initialState.visualize);
    removeFile();
    store.dispatch(updateVisualize(initialState.visualize));
    store.dispatch(updateVisualizeResults(initialState.visualizeResults));
  }

  //   Uploads inputFile and returns a projectID
  async function uploadFile() {
    store.dispatch(
      updateVisualize({
        loading: {
          active: true,
          content: 'Uploading file...',
          showIndicator: true,
        },
      })
    );
    const data = new FormData();
    data.append('file', inputFile);

    let response = await fetch(`${root}visualize/upload`, {
      method: 'POST',
      body: data,
    });

    if (!response.ok) {
      // add error handling
      const { msg, error } = await response.json();
      const message = `<div>
          <p>${msg}</p>
         ${error ? `<p>${error}</p>` : ''} 
        </div>`;
      store.dispatch(updateError({ visible: true, message: message }));
    } else {
      return await response.json();
    }
  }

  function removeFile() {
    setInput(new File([], ''));
  }

  return (
    
    <Form className="mb-2">
    
      <Group controlId="fileType">
        <Label>Choose File Type</Label>
        <Control
          as="select"
          value={inputFormat}
          onChange={(e) => {
            store.dispatch(updateVisualize({ inputFormat: e.target.value }));
          }}
          disabled={disableParameters}
        >
          <option value="vcf">VCF</option>
          <option value="csv">CSV</option>
          <option value="tsv">TSV</option>
          <option value="catalog_tsv">CATALOG TSV</option>
          <option value="catalog_csv">CATALOG CSV</option>
        </Control>
      </Group>
      <Group controlId="fileUpload">
        <section>
          <div {...getRootProps({ className: 'dropzone' })}>
            <input
              {...getInputProps()}
              disabled={inputFile.size || disableParameters}
            />
            {inputFile.size || loading.active ? (
              <button
                id="removeFile"
                className="d-flex w-100 faButton"
                onClick={() => removeFile()}
                disabled={disableParameters}
              >
               
                <span id="uploadedFile">{loading.active ? storeFile : inputFile.name}</span>
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
            store.dispatch(updateVisualize({ selectedGenome: e.target.value }))
          }
          disabled={disableParameters}
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
              store.dispatch(
                updateVisualize({ experimentalStrategy: e.target.value })
              )
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
              store.dispatch(
                updateVisualize({ experimentalStrategy: e.target.value })
              )
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
              store.dispatch(updateVisualize({ mutationSplit: e.target.value }))
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
              store.dispatch(updateVisualize({ mutationSplit: e.target.value }))
            }
            disabled={disableParameters}
          />
          <Check.Label className="font-weight-normal">True</Check.Label>
        </Check>
      </Group>
      <Group className="d-flex">
        <Label className="mr-auto">Multiple Samples</Label>
        <Check inline id="multipleFalse">
          <Check.Input
            type="radio"
            // value={false}
            checked={!isMultiple}
            onChange={() =>
              store.dispatch(updateVisualize({ isMultiple: false }))
            }
            disabled={disableParameters}
          />
          <Check.Label htmlFor="multipleFalse" className="font-weight-normal">
            False
          </Check.Label>
        </Check>
        <Check inline id="multipleTrue">
          <Check.Input
            type="radio"
            // value={true}
            checked={isMultiple}
            onChange={() =>
              store.dispatch(updateVisualize({ isMultiple: true }))
            }
            disabled={disableParameters}
          />
          <Check.Label htmlFor="multipleTrue" className="font-weight-normal">
            True
          </Check.Label>
        </Check>
      </Group>
      <Group className="d-flex">
        <Label className="mr-auto">Collapse Sample</Label>
        <Check inline id="radioFalse">
          <Check.Input
            disabled={!isMultiple || disableParameters}
            type="radio"
            value="False"
            checked={collapseSample == 'False'}
            onChange={(e) =>
              store.dispatch(
                updateVisualize({ collapseSample: e.target.value })
              )
            }
          />
          <Check.Label className="font-weight-normal">False</Check.Label>
        </Check>
        <Check inline id="radioTrue">
          <Check.Input
            disabled={!isMultiple || disableParameters}
            type="radio"
            value="True"
            checked={collapseSample == 'True'}
            onChange={(e) =>
              store.dispatch(
                updateVisualize({ collapseSample: e.target.value })
              )
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
            store.dispatch(updateVisualize({ mutationFilter: e.target.value }))
          }
          disabled={disableParameters}
        ></Control>
      </Group>
      <hr />
      <Group controlId="email">
        <LoadingOverlay active={true} content={'Work in progress...'} />
        <div>
          <Check
            disabled
            inline
            id="toggleQueue"
            type="checkbox"
            label="Submit this job to a Queue"
            checked={queueMode == true}
            onChange={(_) => {
              store.dispatch(updateVisualize({ queueMode: !queueMode }));
            }}
          />
        </div>
        <div>
          <Control
            placeholder="Enter Email"
            size="sm"
            value={email}
            onChange={(e) =>
              store.dispatch(updateVisualize({ email: e.target.value }))
            }
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
      <div className="row">
        <div className="col-sm-6">
          <Button
            disabled={loading.active}
            className="w-100"
            variant="secondary"
            onClick={(e) => handleReset()}
          >
            Reset
          </Button>
        </div>
        <div className="col-sm-6">
          <Button
            disabled={!inputFile.size || disableParameters || loading.active}
            className="w-100"
            variant="primary"
            type="button"
            onClick={(e) => handleSubmit(e)}
          >
            Submit
          </Button>
        </div>
      </div>
    </Form>
  );
}
