import React from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';

export default function TMB({ calculateTMB }) {
  const exploration = useSelector((state) => state.exploration);
  const { plotPath, debugR, err, loading } = exploration.tmb;
  const { projectID, source } = exploration.exposure;

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
        {plotPath && (
          <>
            <hr />
            <Plot
              className="p-3"
              title="Tumor Mutational Burden"
              downloadName={plotPath.split('/').slice(-1)[0]}
              plotPath={`api/results/${projectID}${plotPath}`}
            />
          </>
        )}
      </div>

      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
