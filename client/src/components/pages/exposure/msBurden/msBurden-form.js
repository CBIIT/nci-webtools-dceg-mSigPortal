import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import Select from '../../../controls/select/selectForm';
import Description from '../../../controls/description/description';
import { useSelector, useDispatch } from 'react-redux';
import { actions as exposureActions } from '../../../../services/store/exposure';

import { NavHashLink } from 'react-router-hash-link';
import { customStyles } from '../../../controls/custom/customFormStyle';
const actions = { ...exposureActions };

export default function MsBurdenForm() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.exposure);

  const mergeMsBurden = (state) =>
    dispatch(actions.mergeExposure({ msBurden: state }));

  const { signatureNames } = store.main;
  const { signatureName } = store.msBurden;

  const signatureNameOptions = signatureNames.length
    ? [...new Set(signatureNames.map((d) => d))]
        .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
        .map((e) => ({
          label: e,
          value: e,
        }))
    : [];

  const { control, setValue, watch } = useForm({
    defaultValues: signatureName,
  });

  // set inital
  useEffect(() => {
    if (!signatureName && signatureNameOptions.length) {
      mergeMsBurden({ signatureName: signatureNameOptions[0] });
    }
  }, [signatureNameOptions]);

  return (
    <div>
      <hr />
      <Form className="p-3">
        <Row>
          <Col lg="auto">
            <Select
              name="signatureName"
              label="Signature Name"
              value={signatureName}
              control={control}
              options={signatureNameOptions}
              onChange={(name) => mergeMsBurden({ signatureName: name })}
              styles={customStyles}
              //onChange={handleSignatureName}
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
