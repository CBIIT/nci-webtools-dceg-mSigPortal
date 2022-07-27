import React, { useEffect, useState } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import { actions as exposureActions } from '../../../../services/store/exposure';
import { useSelector, useDispatch } from 'react-redux';

const actions = { ...exposureActions };
const { Group, Label, Check, Control } = Form;
export default function MsPrevalenceForm() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.exposure);
  const mergeMsPrevalence = (state) =>
    dispatch(actions.mergeExposure({ msPrevalence: state }));

  const { mutation } = store.main;
  console.log(mutation);

  const [invalidMin, setMin] = useState(false);

  return (
    <div>
      <hr />
      <Form className="p-3">
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
                type="num"
                name="minnum"
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
          <Col lg="auto" className="d-flex"></Col>
        </Row>
      </Form>
    </div>
  );
}
