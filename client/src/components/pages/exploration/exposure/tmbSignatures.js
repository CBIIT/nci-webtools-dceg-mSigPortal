import React from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';

export default function TmbSignatures({ calculateTmbSig }) {
  const exploration = useSelector((state) => state.exploration);
  const { plotPath, debugR, err, loading } = exploration.tmbSignatures;
  const { projectID, source } = exploration.exposure;

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
              plotPath={`api/results/${projectID}${plotPath}`}
            />
          </>
        )}
      </div>
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
