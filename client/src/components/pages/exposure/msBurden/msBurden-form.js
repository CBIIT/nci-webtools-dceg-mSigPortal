import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import Select from '../../../controls/select/selectForm';
import Description from '../../../controls/description/description';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../../services/store/visualization';
import {
  defaultProfile2,
  defaultMatrix2,
  defaultFilter2,
} from '../../../../services/utils';
import { NavHashLink } from 'react-router-hash-link';

const actions = { ...visualizationActions };

export default function TreeLeafForm() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.exposure);

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
              value={signatureNames}
              options={signatureNameOptions}
              control={control}
            />
          </Col>
          <Col>
            <Button>Recalculate</Button>
          </Col>
        </Row>
      </Form>
    </div>
  );
}
