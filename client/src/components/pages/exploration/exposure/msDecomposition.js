import React from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';

export default function MsDecomposition({ calculateDecomposition }) {
  const dispatch = useDispatch();
  const exploration = useSelector((state) => state.exploration);
  const {
    plotPath,
    txtPath,
    debugR,
    err,
    loading,
  } = exploration.msDecomposition;
  const { projectID, source } = exploration.exposure;

  return (
    <div>
      <div className="p-3">
        <p>
          MS Decomposition: Evaluating the Performance of Mutational Signature
          Decomposition
        </p>
        <p className="m-0">
          This distribution plot below illustrates mutational signature
          decomposition distribution in selected cancer type (by selecting
          Cancer Type Only on the left panel) or across different cancer type.
          Five different methods are used to measure the similarities of
          original mutational profile matrix and reconstructed mutational
          profile matrix across all the samples, including cosine similarity,
          100-L1_Norm_%, 100-L2_Norm_%, KL_Divergence, and Pearson Correlation.
          Click here for the detail of these evaluation method.
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
              onClick={calculateDecomposition}
            >
              Calculate
            </Button>
          </Col>
        </Row>
      </Form>
      <div id="decompositionPlot">
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
              title="Evaluating the Performance of Mutational Signature Decomposition"
              downloadName={plotPath.split('/').slice(-1)[0]}
              plotPath={`api/results/${projectID}${plotPath}`}
              txtPath={projectID + txtPath}
            />
          </>
        )}
      </div>
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
