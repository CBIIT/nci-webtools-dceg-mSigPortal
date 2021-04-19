import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { actions as explorationActions } from '../../../../services/store/exploration';
import { actions as modalActions } from '../../../../services/store/modal';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';

const actions = { ...explorationActions, ...modalActions };

export default function TmbSignatures({ calculateTmbSig }) {
  const dispatch = useDispatch();
  const exploration = useSelector((state) => state.exploration);
  const { plotPath, plotURL, debugR, err, loading } = exploration.tmbSignatures;
  const { projectID, source } = exploration.exposure;
  const mergeExploration = (state) =>
    dispatch(actions.mergeExploration({ exploration: state }));
  const mergeTmbSignatures = (state) =>
    dispatch(actions.mergeExploration({ tmbSignatures: state }));
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
      <div className="p-3">
        <p>TMB Signatures: Tumor Mutational Burden Separated by Signatures</p>
        <p className="m-0">
          The bar plot below illustrates the level of tumor mutational burden
          across different signatures from the selected signature set in the
          selected cancer type. Across the top of the plot are the signatures
          (in the selected cancer type) that were exhibited by at least one
          sample. On the y-axis is the number of mutations per Megabase (log10),
          and the x-axis denotes sample numbers. The green number is the total
          number of samples with the selected cancer type (and therefore
          evaluated for each signature). The blue number is the number of
          samples with the selected cancer type that detected the mutational
          signature from the reference signature set.
        </p>
      </div>
      <hr />
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
