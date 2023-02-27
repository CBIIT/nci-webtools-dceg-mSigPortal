import React, { useEffect, useState } from 'react';
import { Form, Row, Col } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import Select from '../../../controls/select/selectForm';
import { useMsBurdenOptionsQuery } from './apiSlice';

export default function MsBurdenForm({ state, form, setForm }) {
  const [params, setParams] = useState('');
  const { data: signatureNameOptions } = useMsBurdenOptionsQuery(params, {
    skip: !params,
  });

  const { control } = useForm({ defaultValues: form });

  // query signature name options
  useEffect(() => {
    const { study, strategy, signatureSetName, cancer, useAllCancer } = state;
    if (study) {
      setParams({
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        ...(!useAllCancer && { cancer: cancer.value }),
      });
    } else if (state.id) {
      setParams({ userId: state.id });
    }
  }, [state]);

  // set inital
  useEffect(() => {
    if (signatureNameOptions) {
      setForm({ signatureName: signatureNameOptions[0] });
    }
  }, [signatureNameOptions]);

  return (
    <Form className="p-3">
      <Row>
        <Col lg="auto">
          <Select
            name="signatureName"
            label="Signature Name"
            value={form.signatureName}
            disabled={!signatureNameOptions}
            control={control}
            options={signatureNameOptions}
            onChange={(name) => setForm({ signatureName: name })}
          />
        </Col>
      </Row>
    </Form>
  );
}
