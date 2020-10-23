import React, { useState, useEffect } from 'react';
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

export default function Tumor({ submitR }) {
  const rootURL = window.location.pathname;
  const {
    study,
    strategy,
    refSignatureSet,
    loading: mainLoading,
  } = useSelector((state) => state.expExposure);
  const {
    cancer,
    cancerOptions,
    varDataPath,
    plotPath,
    plotURL,
    txtPath,
    debugR,
    err,
    loading,
  } = useSelector((state) => state.expLandscape);

  const [vdFile, setFile] = useState(new File([], ''));

  useEffect(() => {
    if (plotPath) setRPlot(plotPath);
  }, [plotPath]);

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
        const response = await fetch(`${rootURL}results/${plotPath}`);
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

  async function handleSubmit() {
    //   dispatchExpLandscape({ loading: true });
    //   const { projectID } = await upload();
    //   if (projectID) {
    //     calculateR('cosineSimilarity', {
    //       projectID: projectID,
    //       study: study,
    //       cancer: cancer,
    //       strategy: strategy,
    //       refSignatureSet: refSignatureSet,
    //     });
    //   }
    //   dispatchExpLandscape({ loading: false });
    // }
    // async function upload() {
    //   try {
    //     const data = new FormData();
    //     data.append('inputFile', vdFile);
    //     let response = await fetch(`${rootURL}upload`, {
    //       method: 'POST',
    //       body: data,
    //     });
    //     if (!response.ok) {
    //       const { msg, error } = await response.json();
    //       const message = `<div>
    //         <p>${msg}</p>
    //        ${error ? `<p>${error}</p>` : ''}
    //       </div>`;
    //       dispatchError(message);
    //     } else {
    //       return await response.json();
    //     }
    //   } catch (err) {
    //     dispatchError(err);
    //   }
  }

  return (
    <div>
      <Form>
        <LoadingOverlay active={loading || mainLoading} />
        <div>
          <Row className="justify-content-center">
            <Col sm="2">
              <Select
                id="landscapeType"
                label="Cancer Type"
                value={cancer}
                options={cancerOptions}
                onChange={(cancer) =>
                  dispatchExpLandscape({
                    cancer: cancer,
                  })
                }
              />
            </Col>
            <Col sm="4">
              <Label>Upload Variable Data</Label>
              <Form.File
                id="variableData"
                label={vdFile.name || 'Upload here'}
                // accept=''
                onChange={(e) => setFile(e.target.files[0])}
                custom
              />
            </Col>
            <Col sm="5" />
            <Col sm="1" className="m-auto">
              <Button variant="primary" onClick={() => handleSubmit()}>
                Calculate
              </Button>
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
