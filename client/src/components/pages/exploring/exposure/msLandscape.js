import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { actions as exploringActions } from '../../../../services/store/exploring';
import { actions as modalActions } from '../../../../services/store/modal';
import Plot from '../../../controls/plot/plot';
import Select from '../../../controls/select/select';
import Debug from '../../../controls/debug/debug';

const actions = { ...exploringActions, ...modalActions };
const { Label, Group } = Form;

export default function MsLandscape({ calculateLandscape, handleVariable }) {
  const dispatch = useDispatch();
  const exploring = useSelector((state) => state.exploring);
  const {
    variableFile,
    plotPath,
    plotURL,
    txtPath,
    debugR,
    err,
    loading,
  } = exploring.msLandscape;
  const { projectID, source } = exploring.exposure;
  const mergeExploring = (state) =>
    dispatch(actions.mergeExploring({ exploring: state }));
  const mergeMsLandscape = (state) =>
    dispatch(actions.mergeExploring({ msLandscape: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  useEffect(() => {
    plotPath ? setRPlot(plotPath) : clearPlot();
  }, [plotPath, err, debugR, projectID]);

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
          mergeMsLandscape({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        mergeError(err.message);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      mergeMsLandscape({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      mergeMsLandscape({ plotURL: '' });
    }
  }

  return (
    <div>
      <Form className="p-3">
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
          <Col />
          <Col lg="2" className="d-flex">
            <Button
              disabled={source == 'user' && !projectID}
              className="ml-auto mb-auto"
              variant="primary"
              onClick={calculateLandscape}
            >
              Calculate
            </Button>
          </Col>
        </Row>
      </Form>
      <div id="exposureLandscapePlot">
        {err && (
          <div>
            <hr />
            <p className="p-3 text-danger">{err}</p>
          </div>
        )}
        <div style={{ display: plotURL ? 'block' : 'none' }}>
          <hr />
          <Plot
            className="p-3"
            title="Landscape of Mutational Signature Activity"
            downloadName={plotPath.split('/').slice(-1)[0]}
            plotURL={plotURL}
            txtPath={projectID + txtPath}
          />
        </div>
      </div>
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
