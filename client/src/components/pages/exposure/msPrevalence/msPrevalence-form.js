import React, { useEffect, useState } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import { actions as exposureActions } from '../../../../services/store/exposure';
import Description from '../../../controls/description/description';
import { useSelector, useDispatch } from 'react-redux';
import { NavHashLink } from 'react-router-hash-link';

const actions = { ...exposureActions };
const { Group, Label, Check, Control } = Form;
export default function MsPrevalenceForm() {
  const dispatch = useDispatch();
  const exposure = useSelector((state) => state.exposure);
  const { mutation, plotPath, debugR, err, loading } = exposure.msPrevalence;
  const { projectID, source } = exposure.main;
  const mergeMsPrevalence = (state) =>
    dispatch(actions.mergeExposure({ msPrevalence: state }));

  const [invalidMin, setMin] = useState(false);

  return (
    <div>
      <hr />
      <Form noValidate className="p-3">
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
              onClick={() => {}}
              disabled={source == 'user' && !projectID}
            >
              Recalculate
            </Button>
          </Col>
        </Row>
      </Form>
    </div>
  );
}
