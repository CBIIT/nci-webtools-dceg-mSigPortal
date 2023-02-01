import React, { useEffect, useState } from 'react';
import { Form, Row, Col } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import Select from '../../../controls/select/selectForm';
import { useSelector, useDispatch } from 'react-redux';
import { actions as exposureActions } from '../../../../services/store/exposure';
import { useMsBurdenOptionsQuery } from './apiSlice';

const actions = { ...exposureActions };

export default function MsBurdenForm() {
  const dispatch = useDispatch();
  const mergeMsBurden = (state) =>
    dispatch(actions.mergeExposure({ msBurden: state }));

  const { publicForm, main, msBurden } = useSelector((state) => state.exposure);

  const [params, setParams] = useState('');
  const { data: signatureNameOptions } = useMsBurdenOptionsQuery(params, {
    skip: !params,
  });

  const { control } = useForm({ defaultValues: msBurden });

  // query signature name options
  useEffect(() => {
    if (publicForm.study) {
      setParams({
        study: publicForm.study.value,
        strategy: publicForm.strategy.value,
        signatureSetName: publicForm.signatureSetName.value,
        ...(!publicForm.useAllCancer && { cancer: publicForm.cancer.value }),
      });
    } else if (main.id) {
      setParams({ userId: main.id });
    }
  }, [publicForm, main]);

  // set inital
  useEffect(() => {
    if (!msBurden.signatureName && signatureNameOptions) {
      mergeMsBurden({ signatureName: signatureNameOptions[0] });
    }
  }, [signatureNameOptions, msBurden]);

  return (
    <div>
      <hr />
      <Form className="p-3">
        <Row>
          <Col lg="auto">
            <Select
              name="signatureName"
              label="Signature Name"
              value={msBurden.signatureName}
              disabled={!signatureNameOptions}
              control={control}
              options={signatureNameOptions}
              onChange={(name) => mergeMsBurden({ signatureName: name })}
              //onChange={handleSignatureName}
            />
          </Col>
        </Row>
      </Form>
    </div>
  );
}
