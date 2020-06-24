import React, { useState, useCallback } from 'react';
import { Form, Button } from 'react-bootstrap';
import { useDropzone } from 'react-dropzone';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faCloudUploadAlt, faMinus } from '@fortawesome/free-solid-svg-icons';
import './visualize.scss';

const { Group, Label, Control, Check, Text } = Form;

export default function UploadForm() {
  const [inputFormat, setInputFormat] = useState('vcf');
  const [inputFile, setInput] = useState(new File([], ''));
  const [selectedGenome, setSelectedGenome] = useState('GRCh37');
  const [experimentalStrategy, setStrategy] = useState('WGS');
  const [collapseSample, setCollapse] = useState('False');
  const [mutationFilter, setFilter] = useState('text');
  const [mutationSplit, setSplit] = useState('False');
  const [queueMode, setQueueMode] = useState(false);
  const [email, setEmail] = useState('');

  const onDrop = useCallback((acceptedFiles) => {
    setInput(acceptedFiles[0]);
  }, []);
  const { getRootProps, getInputProps } = useDropzone({
    onDrop,
    accept: '.csv, .tsv, .vcf, .gz, .zip, .tar, .tar.gz',
  });

  function handleInputSelect(type) {
    switch (type) {
      case 'vcf':
        setInputFormat('vcf');
        break;
      case 'csv':
        setInputFormat('csv');
        break;
      case 'tsv':
        setInputFormat('tsv');
        break;
      case 'catalog_tsv':
        setInputFormat('catalog_tsv');
        break;
      case 'catalog_csv':
        setInputFormat('catalog_csv');
        break;
      default:
        setInputFormat('vcf');
    }
  }

  function handleGenomeSelect(genome) {
    setSelectedGenome(genome);
  }

  function handleStrategyRadio(experimentalStrategy) {
    setStrategy(experimentalStrategy);
  }

  function handleCollapseRadio(bool) {
    setCollapse(bool);
  }

  function handleFilterInput(string) {
    setFilter(string.trim());
  }

  function handleMutationRadio(bool) {
    setSplit(bool);
  }

  function handleEmailInput(string) {
    setEmail(string.trim());
  }

  function submitHandler() {
    const upload = uploadFile();
    upload.then((data) => {
      if (data && data.projectID) {
        let req = {
          projectID: data.projectID,
          inputFile: inputFile,
          inputFormat: inputFormat,
          genomeAssemblyVersion: selectedGenome,
          experimentalStrategy: experimentalStrategy,
          collapseSample: collapseSample,
          mutationFilter: mutationFilter,
          mutationSplit: mutationSplit,
        };
        console.log('fn args', req);
      } else {
        //   add error handling
      }
    });
  }

  //   Uploads inputFile and returns a projectID
  async function uploadFile() {
    const root =
      process.env.NODE_ENV === 'development'
        ? 'http://localhost:8330/'
        : window.location.pathname;
    const data = new FormData();
    data.append('file', inputFile);

    let response = await fetch(`${root}upload`, {
      method: 'POST',
      body: data,
    });

    if (!response.ok) {
      // add error handling
      console.log(`HTTP error! status: ${response.status}`);
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
          onChange={(e) => handleInputSelect(e.target.value)}
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
            <input {...getInputProps()} disabled={inputFile.size} />
            {inputFile.size ? (
              <button
                id="removeFile"
                className="d-flex w-100"
                onClick={() => removeFile()}
              >
                <span id="uploadedFile">{inputFile.name}</span>
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
          onChange={(e) => handleGenomeSelect(e.target.value)}
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
            onChange={(e) => handleStrategyRadio(e.target.value)}
          />
          <Check.Label className="font-weight-normal">WGS</Check.Label>
        </Check>
        <Check inline id="radioWES">
          <Check.Input
            type="radio"
            value="WES"
            checked={experimentalStrategy == 'WES'}
            onChange={(e) => handleStrategyRadio(e.target.value)}
          />
          <Check.Label className="font-weight-normal">WES</Check.Label>
        </Check>
      </Group>
      <Group className="d-flex">
        <Label className="mr-auto">Collapse Sample</Label>
        <Check inline id="radioFalse">
          <Check.Input
            type="radio"
            value="False"
            checked={collapseSample == 'False'}
            onChange={(e) => handleCollapseRadio(e.target.value)}
          />
          <Check.Label className="font-weight-normal">False</Check.Label>
        </Check>
        <Check inline id="radioTrue">
          <Check.Input
            type="radio"
            value="True"
            checked={collapseSample == 'True'}
            onChange={(e) => handleCollapseRadio(e.target.value)}
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
          onChange={(e) => handleFilterInput(e.target.value)}
        ></Control>
        <Text className="text-muted">Default: text</Text>
      </Group>
      <Group className="d-flex">
        <Label className="mr-auto">Mutation Split</Label>
        <Check inline id="radioMutationSplitFalse">
          <Check.Input
            type="radio"
            value="False"
            checked={mutationSplit == 'False'}
            onChange={(e) => handleMutationRadio(e.target.value)}
          />
          <Check.Label className="font-weight-normal">False</Check.Label>
        </Check>
        <Check inline id="radioMutationSplitTrue">
          <Check.Input
            type="radio"
            value="True"
            checked={mutationSplit == 'True'}
            onChange={(e) => handleMutationRadio(e.target.value)}
          />
          <Check.Label className="font-weight-normal">True</Check.Label>
        </Check>
      </Group>
      <Group controlId="email">
        <div>
          <Check
            disabled
            inline
            id="toggleQueue"
            type="checkbox"
            label="Submit this job to a Queue"
            checked={queueMode == true}
            onChange={(_) => {
              setQueueMode(!queueMode);
            }}
          />
        </div>
        <div>
          <Control
            placeholder="Enter Email"
            size="sm"
            value={email}
            onChange={(e) => handleEmailInput(e.target.value)}
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
          <Button className="w-100" variant="secondary">
            Reset
          </Button>
        </div>
        <div className="col-sm-6">
          <Button
            disabled={!inputFile.size}
            className="w-100"
            variant="primary"
            type="submit"
            onClick={() => submitHandler()}
          >
            Submit
          </Button>
        </div>
      </div>
    </Form>
  );
}
