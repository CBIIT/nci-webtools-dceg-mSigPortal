import React, { useState } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { actions as exposureActions } from '../../../services/store/exposure';
import { actions as modalActions } from '../../../services/store/modal';
import Plot from '../../controls/plot/plot';
import Debug from '../../controls/debug/debug';

const actions = { ...exposureActions, ...modalActions };
const { Group, Label, Control } = Form;

export default function MsPrevalence({ calculatePrevalence }) {
  const dispatch = useDispatch();
  const exposure = useSelector((state) => state.exposure);
  const { mutation, plotPath, debugR, err, loading } = exposure.msPrevalence;
  const { projectID, source } = exposure.exposureState;
  const mergeMsPrevalence = (state) =>
    dispatch(actions.mergeExposure({ msPrevalence: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const [invalidMin, setMin] = useState(false);

  return (
    <div>
      <div className="p-3">
        <p>MS Prevalence: Prevalence of Mutational Signature </p>
        <p>
          This page allows you to analyze the prevalence of signatures from the
          selected Study and Reference Signature Set by both sample and
          mutation. For prevalence by samples, Input the “Minimal Number of
          Mutations Assigned to Each Signature” to set the smallest number of
          mutations assigned to each signature within sample, which can have to
          be included in the result.
        </p>
        <p className="m-0">
          The pie chart on the left illustrates the prevalence of each
          mutational signature by mutations. The bar plot on the right
          illustrates the prevalence of each mutational signature by samples.
          The colors represent each of the mutational signatures in both plots.
        </p>
      </div>
      <hr />
      <Form noValidate className="p-3">
        <LoadingOverlay active={loading} />
        <Row>
          <Col lg="5">
            <Group
              controlId="prevalenceMutations"
              title="Minimum Number of Mutations Within Each Signature"
            >
              <Label>
                Minimum Number of Mutations Assigned to Each Signature
              </Label>
              <Control
                value={mutation}
                placeholder="e.g. 100"
                onChange={(e) => {
                  mergeMsPrevalence({
                    mutation: e.target.value,
                  });
                }}
                isInvalid={invalidMin}
              />
              <Form.Control.Feedback type="invalid">
                Enter a numeric minimum value
              </Form.Control.Feedback>
            </Group>
          </Col>
          <Col />
          <Col lg="2" className="d-flex">
            <Button
              className="ml-auto mb-auto"
              variant="primary"
              onClick={() => {
                if (!mutation || isNaN(mutation)) setMin(true);
                else setMin(false);

                if (mutation && !isNaN(mutation)) calculatePrevalence();
              }}
              disabled={source == 'user' && !projectID}
            >
              Calculate
            </Button>
          </Col>
        </Row>
      </Form>
      <div id="exposurePrevalencePlot">
        {err && (
          <div>
            <p className="p-3 text-danger">{err}</p>
          </div>
        )}
        {plotPath && (
          <>
            <hr />
            <Plot
              className="p-3"
              title="Prevalence of Mutational Signature"
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
