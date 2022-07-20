import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import Select from '../../../controls/select/selectForm';
import Description from '../../../controls/description/description';
import { useSelector, useDispatch } from 'react-redux';
import { actions as exposureActions } from '../../../../services/store/exposure';

import { NavHashLink } from 'react-router-hash-link';

const actions = { ...exposureActions };

export default function TreeLeafForm() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.exposure);

  const mergeMsBurden = (state) =>
    dispatch(actions.mergeExposure({ msBurden: state }));

  const { signatureNames } = store.main;

  console.log(signatureNames);
  const { control, setValue, watch } = useForm();

  const signatureNameOptions = signatureNames.length
    ? [...new Set(signatureNames.map((d) => d))]
        .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
        .map((e) => ({
          label: e,
          value: e,
        }))
    : [];

  return (
    <div>
      <hr />
      <Form className="p-3">
        <Row>
          <Col lg="auto">
            <Select
              name="signatureNames"
              label="Signature Name"
              control={control}
              options={signatureNameOptions}
              onChange={(name) => mergeMsBurden({ signatureName: name })}
            />
          </Col>
          {/* <Col lg="auto" className="d-flex">
            <Button
              className="mt-auto mb-2"
              variant="primary"
              //onClick={}
            >
              Recalculate
            </Button>
          </Col> */}
        </Row>
      </Form>
    </div>
  );
}
