import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchExpLandscape,
} from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';

const { Label, Group } = Form;

export default function Landscape({ calculateLandscape, handleVariable }) {
  const {
    variableFile,
    plotPath,
    plotURL,
    txtPath,
    debugR,
    err,
    loading,
  } = useSelector((state) => state.expLandscape);
  const { projectID } = useSelector((state) => state.expExposure);

  useEffect(() => {
    if (plotPath) setRPlot(plotPath);
    else clearPlot();
  }, [plotPath]);

  async function setRPlot(plotPath) {
    if (plotPath) {
      try {
        const response = await fetch(`api/results/${projectID}${plotPath}`);
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

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      dispatchExpLandscape({ plotPath: '', plotURL: '' });
    }
  }

  return (
    <div>
      <Form>
        <LoadingOverlay active={loading} />
        <Row className="">
          <Col lg="4">
            <Group controlId="landscape">
              <Label>Upload Variable Data</Label>
              <Form.File
                id="variableData"
                label={variableFile || 'Upload here (optional)'}
                // accept=''
                onChange={(e) => handleVariable(e.target.files[0])}
                custom
              />
            </Group>
          </Col>
          <Col lg="1" className="d-flex justify-content-end">
            <Button
              className="mt-auto mb-3"
              variant="secondary"
              onClick={() => {
                handleVariable(new File([], ''));
              }}
              disabled={!variableFile}
            >
              Remove
            </Button>
          </Col>
          <Col lg="5" />
          <Col lg="2" className="d-flex justify-content-end">
            <Button
              className="mt-auto mb-3"
              variant="primary"
              onClick={calculateLandscape}
            >
              Calculate
            </Button>
          </Col>
        </Row>
        <div id="exposureLandscapePlot">
          {err && (
            <div>
              <p>
                An error has occured. Please verify your input.
              </p>
              <p>Error: {err}</p>
            </div>
          )}
          <div style={{ display: plotURL ? 'block' : 'none' }}>
            <Plot
              plotName={plotPath.split('/').slice(-1)[0]}
              plotURL={plotURL}
              txtPath={projectID + txtPath}
            />
          </div>
        </div>
      </Form>
      <Debug msg={debugR} />
    </div>
  );
}
