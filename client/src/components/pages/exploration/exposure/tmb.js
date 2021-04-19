import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { actions as explorationActions } from '../../../../services/store/exploration';
import { actions as modalActions } from '../../../../services/store/modal';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';

const actions = { ...explorationActions, ...modalActions };

export default function TMB({ calculateTMB }) {
  const dispatch = useDispatch();
  const exploration = useSelector((state) => state.exploration);
  const { plotPath, plotURL, debugR, err, loading } = exploration.tmb;
  const { projectID, source } = exploration.exposure;
  const mergeExploration = (state) =>
    dispatch(actions.mergeExploration({ exploration: state }));
  const mergeTMB = (state) =>
    dispatch(actions.mergeExploration({ tmb: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  useEffect(() => {
    plotPath ? setRPlot(plotPath) : clearPlot();
  }, [plotPath, err, debugR]);

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
          mergeTMB({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        mergeError(err.message);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      mergeTMB({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      mergeTMB({ plotURL: '' });
    }
  }

  return (
    <div>
      <div className="p-3">
        <p> TMB: Tumor Mutational Burden</p>
        <p className="m-0">
          The bar plot below illustrates the level of tumor mutational burden
          (number of mutations per megabase) across different cancer types for
          selected study. Across the top of the plot are the different cancer
          types. The y-axis is the number of mutations per megabase (log10), and
          the x-axis denotes sample numbers. The green number is the number of
          samples for a given cancer type, and the blue number is the number of
          samples that had mutation data for that cancer type.
        </p>
      </div>
      <hr />
      <Form className="p-3">
        <LoadingOverlay active={loading} />
        <Row>
          <Col>Select parameters from the left side panel.</Col>
          <Col lg="2" className="d-flex">
            <Button
              disabled={source == 'user' && !projectID}
              className="ml-auto mb-auto"
              variant="primary"
              onClick={calculateTMB}
            >
              Calculate
            </Button>
          </Col>
        </Row>
      </Form>
      <div id="tmbPlot">
        {err && (
          <div>
            <hr />
            <p className="p-3 text-danger">{err}</p>
          </div>
        )}
        {plotURL && (
          <>
            <hr />
            <Plot
              className="p-3"
              title="Tumor Mutational Burden"
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
