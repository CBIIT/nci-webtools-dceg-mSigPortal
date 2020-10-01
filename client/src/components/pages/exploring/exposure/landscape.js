import React, { useState } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchExpLandscape,
} from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';
import Select from '../../../controls/select/select';

const { Group, Label, Control } = Form;

export default function Tumor({ submitR, downloadResults }) {
  const rootURL = window.location.pathname;
  const {
    study,
    studyOptions,
    cancer,
    cancerOptions,
    strategy,
    strategyOptions,
    refSignatureSet,
    refSignatureSetOptions,
    genomeSize,
    plotPath,
    plotURL,
    txtPath,
    debugR,
    err,
    displayDebug,
    loading,
  } = useSelector((state) => state.expLandscape);
  const { displayTab, publicDataOptions } = useSelector(
    (state) => state.exploring
  );

  const [vdFile, setFile] = useState(new File([], ''));

  async function calculateR(fn, args) {
    console.log(fn);
    dispatchExpLandscape({
      loading: true,
      err: false,
      debugR: '',
    });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();

        dispatchExpLandscape({
          loading: false,
          debugR: err,
        });
      } else {
        const { debugR, output } = await response.json();

        dispatchExpLandscape({
          debugR: debugR,
          loading: false,
          plotPath: output.plotPath,
          txtPath: output.txtPath,
        });
        setRPlot(output.plotPath);
      }
    } catch (err) {
      dispatchError(err);
      dispatchExpLandscape({ loading: false });
    }
  }

  async function setRPlot(plotPath) {
    if (plotPath) {
      try {
        const response = await fetch(`${rootURL}getSVG`, {
          method: 'POST',
          headers: {
            Accept: 'image/svg',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({ path: plotPath }),
        });
        if (!response.ok) {
          // console.log(await response.json());
        } else {
          const pic = await response.blob();
          const objectURL = URL.createObjectURL(pic);

          if (plotURL) URL.revokeObjectURL(plotURL);
          dispatchExpLandscape({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchExpLandscape({ err: true, plotURL: '' });
    }
  }

  function handleStudy(study) {
    const cancerOptions = [
      ...new Set(
        publicDataOptions
          .filter((data) => data.Study == study)
          .map((data) => data.Cancer_Type)
      ),
    ];
    const strategyOptions = [
      ...new Set(
        publicDataOptions
          .filter(
            (data) =>
              data.Study == study && data.Cancer_Type == cancerOptions[0]
          )
          .map((data) => data.Dataset)
      ),
    ];

    dispatchExpLandscape({
      study: study,
      cancer: cancerOptions[0],
      cancerOptions: cancerOptions,
      strategy: strategyOptions[0],
      strategyOptions: strategyOptions,
    });
  }

  function handleCancer(cancer) {
    const strategyOptions = [
      ...new Set(
        publicDataOptions
          .filter((data) => data.Study == study && data.Cancer_Type == cancer)
          .map((data) => data.Dataset)
      ),
    ];

    dispatchExpLandscape({
      cancer: cancer,
      strategy: strategyOptions[0],
      strategyOptions: strategyOptions,
    });
  }

  async function handleSubmit() {
    dispatchExpLandscape({ loading: true });

    const { projectID } = await upload();
    if (projectID) {
      calculateR('cosineSimilarity', {
        projectID: projectID,
        study: study,
        cancer: cancer,
        strategy: strategy,
        refSignatureSet: refSignatureSet,
        genomeSize: parseFloat(genomeSize),
      });
    }
    dispatchExpLandscape({ loading: false });
  }

  async function upload() {
    try {
      const data = new FormData();
      data.append('inputFile', vdFile);
      let response = await fetch(`${rootURL}upload`, {
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
    }
  }

  return (
    <div>
      <Form>
        <LoadingOverlay active={loading} />
        <div>
          <Row className="justify-content-center">
            <Col sm="2">
              <Select
                id="landscapeStudy"
                label="Study"
                value={study}
                options={studyOptions}
                onChange={(study) => handleStudy(study)}
              />
            </Col>
            <Col sm="2">
              <Select
                id="landscapeType"
                label="Cancer Type"
                value={cancer}
                options={cancerOptions}
                onChange={(cancer) => handleCancer(cancer)}
              />
            </Col>
            <Col sm="2">
              <Select
                id="landscapeStrategy"
                label="Experimental Strategy"
                value={strategy}
                options={strategyOptions}
                onChange={(strategy) =>
                  dispatchExpLandscape({ strategy: strategy })
                }
              />
            </Col>
            <Col sm="3">
              <Select
                id="landscapeRefSet"
                label="Reference Signature Set"
                value={refSignatureSet}
                options={refSignatureSetOptions}
                onChange={(set) =>
                  dispatchExpLandscape({ refSignatureSet: set })
                }
              />
            </Col>
            <Col sm="2">
              <Group controlId="landscapeGenomeSize">
                <Label>Genome Size</Label>
                <Control
                  id="landscapeGenomeSize"
                  value={genomeSize}
                  onChange={(e) => {
                    dispatchExpLandscape({
                      genomeSize: e.target.value,
                    });
                  }}
                />
                {/* <Text className="text-muted">(Ex. NCG>NTG)</Text> */}
              </Group>
            </Col>
            <Col sm="1" className="m-auto">
              <Button variant="primary" onClick={() => handleSubmit()}>
                Calculate
              </Button>
            </Col>
          </Row>
          <Row>
            <Col sm="4">
              <Label>Upload Variable Data</Label>
              <Form.File
                id="variableData"
                label={vdFile.name || 'Upload here'}
                onChange={(e) => setFile(e.target.files[0])}
                custom
              />
            </Col>
          </Row>
          <div id="withinPlot">
            <div style={{ display: err ? 'block' : 'none' }}>
              <p>
                An error has occured. Check the debug section for more info.
              </p>
            </div>
            <div style={{ display: plotURL ? 'block' : 'none' }}>
              <Plot
                plotName={plotPath.split('/').slice(-1)[0]}
                plotURL={plotURL}
                txtPath={txtPath}
              />
            </div>
          </div>
        </div>
      </Form>
      <Debug msg={debugR} />
    </div>
  );
}
