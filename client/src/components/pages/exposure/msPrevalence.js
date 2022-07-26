import React, { useState } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import Description from '../../controls/description/description';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { actions as exposureActions } from '../../../services/store/exposure';
import SvgContainer from '../../controls/svgContainer/svgContainer';
import Debug from '../../controls/debug/debug';

const actions = { ...exposureActions };
const { Group, Label, Control } = Form;

export default function MsPrevalence({ calculatePrevalence }) {
  const dispatch = useDispatch();
  const exposure = useSelector((state) => state.exposure);
  const { mutation, plotPath, debugR, err, loading } = exposure.msPrevalence;
  const { projectID, source } = exposure.main;
  const mergeMsPrevalence = (state) =>
    dispatch(actions.mergeExposure({ msPrevalence: state }));

  const [invalidMin, setMin] = useState(false);

  return (
    <div>
      <div className="p-3">
        <b>Prevalence of Mutational Signature </b>
        <Description
          less="The following plot indicates both mutation and sample level prevalence of signatures from the selected Study."
          more="For prevalence by samples, input the [Minimal Number of Mutations Assigned to Each Signature] to set the smallest number of mutations assigned to each signature required for the detection of the mutational signature in each sample."
        />
      </div>
      <hr />
      <Form noValidate className="p-3">
        <LoadingOverlay active={loading} />
        <Row>
          <Col lg="auto">
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
          <Col lg="auto" className="d-flex">
            <Button
              className="mt-auto mb-3"
              variant="primary"
              onClick={() => {
                if (!mutation || isNaN(mutation)) setMin(true);
                else setMin(false);

                if (mutation && !isNaN(mutation)) calculatePrevalence();
              }}
              disabled={source == 'user' && !projectID}
            >
              Recalculate
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
            <SvgContainer
              className="p-3"
              title="Prevalence of Mutational Signatures"
              downloadName={plotPath.split('/').slice(-1)[0]}
              plotPath={`results/${plotPath}`}
            />
            <p className="p-3">
              The pie chart on the left illustrates the prevalence of each
              mutational signature by mutations. The bar plot on the right
              illustrates the prevalence of each mutational signature by
              samples. The colors represent the mutational signatures in both
              plots.
            </p>
          </>
        )}
      </div>
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
