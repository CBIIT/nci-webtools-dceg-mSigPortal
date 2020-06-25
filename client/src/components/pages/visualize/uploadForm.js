import React, { useState, useCallback } from 'react';
import { Form, Button } from 'react-bootstrap';
import { useDropzone } from 'react-dropzone';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faCloudUploadAlt, faMinus } from '@fortawesome/free-solid-svg-icons';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import './visualize.scss';

const { Group, Label, Control, Check, Text } = Form;

const root =
  process.env.NODE_ENV === 'development'
    ? 'http://localhost:8330/'
    : window.location.pathname;

export default function UploadForm({ setPlots, setOpenSidebar }) {
  const [inputFormat, setInputFormat] = useState('vcf');
  const [inputFile, setInput] = useState(new File([], ''));
  const [selectedGenome, setSelectedGenome] = useState('GRCh37');
  const [experimentalStrategy, setStrategy] = useState('WGS');
  const [mutationSplit, setSplit] = useState('False');
  const [isMultiple, setMultiple] = useState(false);
  const [collapseSample, setCollapse] = useState('False');
  const [mutationFilter, setFilter] = useState('');
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

  function handleMutationRadio(bool) {
    setSplit(bool);
  }

  function handleMultiple(bool) {
    setMultiple(bool);
  }

  function handleCollapseRadio(bool) {
    setCollapse(bool);
  }

  function handleFilterInput(string) {
    setFilter(string.trim());
  }

  function handleEmailInput(string) {
    setEmail(string.trim());
  }

  async function handleSubmit(e) {
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
      if (mutationFilter.length)
        args['mutationFilter'] = ['-F', mutationFilter];
      if (inputFormat == 'vcf') args['mutationSplit'] = ['-s', mutationSplit];
      if (isMultiple) args['collapseSample'] = ['-c', collapseSample];

      const response = await fetch(`${root}api/visualize`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(args),
      });
      const result = await response.json();
      console.log(result);
      setPlots(result);
      setOpenSidebar(false);
    } else {
      //   add error handling
    }
    e.preventDefault();
  }

  //   Uploads inputFile and returns a projectID
  async function uploadFile() {
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
        <Label className="mr-auto">Mutation Split</Label>
        <Check inline id="radioMutationSplitFalse">
          <Check.Input
            disabled={inputFormat != 'vcf'}
            type="radio"
            value="False"
            checked={mutationSplit == 'False'}
            onChange={(e) => handleMutationRadio(e.target.value)}
          />
          <Check.Label className="font-weight-normal">False</Check.Label>
        </Check>
        <Check inline id="radioMutationSplitTrue">
          <Check.Input
            disabled={inputFormat != 'vcf'}
            type="radio"
            value="True"
            checked={mutationSplit == 'True'}
            onChange={(e) => handleMutationRadio(e.target.value)}
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
            onChange={() => handleMultiple(!isMultiple)}
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
            onChange={() => handleMultiple(!isMultiple)}
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
            disabled={!isMultiple}
            type="radio"
            value="False"
            checked={collapseSample == 'False'}
            onChange={(e) => handleCollapseRadio(e.target.value)}
          />
          <Check.Label className="font-weight-normal">False</Check.Label>
        </Check>
        <Check inline id="radioTrue">
          <Check.Input
            disabled={!isMultiple}
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
            onClick={(e) => handleSubmit(e)}
          >
            Submit
          </Button>
        </div>
      </div>
    </Form>
  );
}
