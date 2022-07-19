import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
// { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Select from '../../../controls/select/selectForm';
import Description from '../../../controls/description/description';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../../services/store/visualization';
//import CustomSelect from '../../controls/select/select-old';

import { NavHashLink } from 'react-router-hash-link';

const actions = { ...visualizationActions };

export default function MsBurden() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.exposure);

  const { signatureName } = store.main;
  console.log(signatureName);

  //   const mergeMsBurden = (state) =>
  //     dispatch(actions.mergeExposure({ msBurden: state }));

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
              //   options={
              //     source == 'public' ? signatureNameOptions : userNameOptions
              //   }
              //onChange={(name) => mergeMsBurden({ signatureName: name })}
            />
          </Col>
          <Col lg="auto" className="d-flex">
            <Button
              className="mt-auto mb-2"
              variant="primary"
              //onClick={calculateBurden}
            >
              Recalculate
            </Button>
          </Col>
        </Row>
      </Form>
    </div>
  );
}
