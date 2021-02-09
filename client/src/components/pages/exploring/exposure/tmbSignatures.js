import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { actions as exploringActions } from '../../../../services/store/exploring';
import { actions as modalActions } from '../../../../services/store/modal';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';

const actions = { ...exploringActions, ...modalActions };

export default function TmbSignatures({ calculateTmbSig }) {
  const dispatch = useDispatch();
  const exploring = useSelector((state) => state.exploring);
  const { plotPath, plotURL, debugR, err, loading } = exploring.tmbSignatures;
  const { projectID, source } = exploring.exposure;
  const mergeExploring = (state) =>
    dispatch(actions.mergeExploring({ exploring: state }));
  const mergeTmbSignatures = (state) =>
    dispatch(actions.mergeExploring({ tmbSignatures: state }));
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
          mergeTmbSignatures({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        mergeError(err.message);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      mergeTmbSignatures({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      mergeTmbSignatures({ plotURL: '' });
    }
  }

  return (
    <div>
      <Form className="p-3">
        <LoadingOverlay active={loading} />
        <Row>
          <Col>Select parameters from the left side panel.</Col>
          <Col />
          <Col lg="2" className="d-flex">
            <Button
              disabled={source == 'user' && !projectID}
              className="ml-auto mb-auto"
              variant="primary"
              onClick={calculateTmbSig}
            >
              Calculate
            </Button>
          </Col>
        </Row>
      </Form>
      <div id="tmbSigPlot">
        {err && (
          <div>
            <hr />
            <p className="p-3 text-danger">{err}</p>
          </div>
        )}
        {plotPath && (
          <>
            <hr />
            <Plot
              className="p-3"
              title="Tumor Mutational Burden Separated by Signatures"
              downloadName={plotPath.split('/').slice(-1)[0]}
              plotURL={plotURL}
            />
          </>
        )}
      </div>
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
