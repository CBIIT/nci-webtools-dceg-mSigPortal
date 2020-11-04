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

const { Label } = Form;

export default function Landscape({ calculateLandscape }) {
  const rootURL = window.location.pathname;
  const { loading: mainLoading } = useSelector((state) => state.expExposure);
  const { plotPath, plotURL, txtPath, debugR, err, loading } = useSelector(
    (state) => state.expLandscape
  );

  const [vdFile, setFile] = useState(new File([], ''));

  useEffect(() => {
    if (plotPath) setRPlot(plotPath);
  }, [plotPath]);

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

  async function handleUpload() {
    if (vdFile.size) {
      dispatchExpLandscape({ loading: true });

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
          const { filePath } = await response.json();
          dispatchExpLandscape({ varDataPath: filePath });
        }
      } catch (err) {
        dispatchError(err);
      }

      dispatchExpLandscape({ loading: false });
    }
  }

  return (
    <div>
      <Form>
        <LoadingOverlay active={loading || mainLoading} />
        <div className="px-4">
          <Row className="justify-content-center">
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
            <Col sm="1" className="my-auto">
              <Button
                variant="secondary"
                onClick={handleUpload}
                disabled={!vdFile.size}
              >
                Upload
              </Button>
            </Col>
            <Col sm="6" />
            <Col sm="1" className="m-auto">
              <Button variant="primary" onClick={calculateLandscape}>
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
