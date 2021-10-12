import React from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { actions as exposureActions } from '../../../services/store/exposure';
import { actions as modalActions } from '../../../services/store/modal';
import Plot from '../../controls/plot/plot';
import Select from '../../controls/select/select';
import Debug from '../../controls/debug/debug';

const actions = { ...exposureActions, ...modalActions };

export default function MsBurden({ calculateBurden }) {
  const dispatch = useDispatch();
  const exposure = useSelector((state) => state.exposure);

  const { signatureName, plotPath, debugR, err, loading } = exposure.msBurden;
  const { projectID, signatureNameOptions, userNameOptions, source } =
    exposure.exposureState;

  const mergeMsBurden = (state) =>
    dispatch(actions.mergeExposure({ msBurden: state }));

  return (
    <div>
      <div className="p-3">
        <b>Mutational Signature Burden Across Cancer Types</b>
        <p className="m-0">
          The bar plot below illustrates the level of mutational signature
          burden across different cancer types with regard to a specific
          signature. The signature selected is one from the Reference Signature
          Set selected in the left panel. Across the top of the plot are the
          cancer types. On the y-axis is the number of mutations per Megabase
          (log10), and the x-axis denotes sample numbers. The green number is
          the number of samples with the cancer type (across the top of the
          plot), and the blue number is the number of samples in that cancer
          type that detected the signature selected.
        </p>
      </div>
      <hr />
      <Form className="p-3">
        <LoadingOverlay active={loading} />
        <Row>
          <Col lg="auto">
            <Select
              className="mb-2"
              disabled={source == 'user' && !userNameOptions.length}
              id="acrossSignatureName"
              label="Signature Name"
              value={signatureName}
              options={
                source == 'public' ? signatureNameOptions : userNameOptions
              }
              onChange={(name) => mergeMsBurden({ signatureName: name })}
            />
          </Col>
          <Col lg="auto" className="d-flex">
            <Button
              disabled={!signatureName || (source == 'user' && !projectID)}
              className="mt-auto mb-2"
              variant="primary"
              onClick={calculateBurden}
            >
              Recalculate
            </Button>
          </Col>
        </Row>
      </Form>
      <div id="exposureAcrossPlot">
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
              title="Mutational Signature Burden Across Cancer Types"
              downloadName={plotPath.split('/').slice(-1)[0]}
              plotPath={`api/results/${plotPath}`}
            />
          </>
        )}
        {/* <Debug msg={debugR} /> */}
      </div>
    </div>
  );
}
