import React, { useState } from 'react';
import { Form, Button } from 'react-bootstrap';
import { useDropzone } from 'react-dropzone';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faFileUpload } from '@fortawesome/free-solid-svg-icons';
import './visualize.scss';

const { Group, Label, Control, Check, Text, FormCheck } = Form;

export default function UploadForm() {
  const [inputFormat, setInputFormat] = useState('vcf');
  const [selectedFileType, setSelectedFile] = useState('vcf');
  const [projectID, setProjectID] = useState('');
  const [selectedGenome, setSelectedGenome] = useState('GRCh37');
  const [experimentalStrategy, setStrategy] = useState('WGS');
  const [collapseSample, setCollapse] = useState('False');
  const [mutationFilter, setFilter] = useState('text');
  const [mutationSplit, setSplit] = useState('False');

  const { acceptedFiles, getRootProps, getInputProps } = useDropzone({
    accept: '.csv,.tsv.,.vcf,.gz,.zip,.tar,.tar.gz',
  });

  const files = acceptedFiles.map((file) => (
    <li key={file.path}>
      {file.path} - {file.size} bytes
    </li>
  ));

  function handleInputSelect(type) {
    switch (type) {
      case 'vcf':
        setInputFormat('vcf');
        setSelectedFile('vcf');
        break;
      case 'csv':
        setInputFormat('csv');
        setSelectedFile('csv');
        break;
      case 'tsv':
        setInputFormat('tsv');
        setSelectedFile('tsv');
        break;
      case 'cat_tsv':
        setInputFormat('catalog_tsv');
        setSelectedFile('catalog_tsv');
        break;
      case 'cat_csv':
        setInputFormat('catalog_csv');
        setSelectedFile('catalog_csv');
        break;
      default:
        setInputFormat('vcf');
        setSelectedFile('vcf');
    }
  }

  function handleProjectInput(id) {
    setProjectID(id.trim());
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

  return (
    <Form className="mb-2">
      <Group controlId="fileType">
        <Label className="required">Choose File Type</Label>
        <Control
          as="select"
          value={selectedFileType}
          onChange={(e) => handleInputSelect(e.target.value)}
        >
          <option value="vcf">VCF</option>
          <option value="csv">CSV</option>
          <option value="tsv">TSV</option>
          <option value="cat_tsv">CATALOG TSV</option>
          <option value="cat_csv">CATALOG CSV</option>
        </Control>
      </Group>
      <Group controlId="fileUpload">
        <section>
          <div {...getRootProps({ className: 'dropzone' })}>
            <input {...getInputProps()} />
            <p>Drop files here or click to upload.</p>
            <FontAwesomeIcon icon={faFileUpload} size="4x" />
          </div>
          <aside>
            <ul>{files}</ul>
          </aside>
        </section>
      </Group>
      <Group controlId="projectID">
        <Label className="required">Project ID</Label>
        <Control
          type="text"
          size="sm"
          placeholder="Enter a Project ID"
          value={projectID}
          onChange={(e) => handleProjectInput(e.target.value)}
        ></Control>
        <Text className="text-muted">Input Rules: No Whitespace, ... TBA</Text>
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
      <Group controlId="projectID">
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
      <div className="d-flex">
        <Button className="ml-auto" variant="primary" type="submit">
          Submit
        </Button>
      </div>
    </Form>
  );
}
