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
  const { projectID, source } = store.main;
  const mergeMsPrevalence = (state) =>
    dispatch(actions.mergeExposure({ msPrevalence: state }));

  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const { mutation } = store.main;
  console.log(mutation);

  const [invalidMin, setMin] = useState(false);

  async function calculatePrevalence() {
    try {
      if (source == 'user' && !projectID) {
        mergeError('Missing Required Files');
      } else {
        mergeMsPrevalence({
          loading: true,
          err: false,
          plotPath: '',
        });

        //await handleCalculate('prevalence');

        mergeMsPrevalence({ loading: false });
      }
    } catch (error) {
      mergeError(error.message);
    }
  }

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
          <Col lg="auto" className="d-flex">
            <Button
              className="mt-auto mb-3"
              variant="primary"
              onClick={() => {
                calculatePrevalence();
              }}
            >
              Recalculate
            </Button>
          </Col>
        </Row>
      </Form>
    </div>
  );
}
