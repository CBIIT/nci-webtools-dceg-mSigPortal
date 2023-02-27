import React, { useEffect, useState } from 'react';
import { Form, Row, Col } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import Select from '../../../controls/select/selectForm';
import { useMsAssociationOptionsQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

const { Group, Check } = Form;

export default function MsAssociationForm({ state, form, mergeForm }) {
  const [signatureOptionParams, setSignatureOptionParams] = useState('');

  const { data, error, isFetching } = useMsAssociationOptionsQuery(
    signatureOptionParams,
    {
      skip: !signatureOptionParams,
    }
  );

  const { both } = form;

  const { study, strategy, signatureSetName, cancer, useAllCancer, id } = state;

  const { control } = useForm({
    defaultValues: form,
  });

  // set inital
  useEffect(() => {
    if (data) {
      mergeForm({
        signatureName1: data[0],
        signatureName2: data[1] || data[0],
      });
    }
  }, [data]);

  // get signature name options
  useEffect(() => {
    if (study) {
      setSignatureOptionParams({
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        ...(!useAllCancer && { cancer: cancer.value }),
      });
    } else if (id) {
      setSignatureOptionParams({ userId: id });
    }
  }, [state]);

  function handleSignatureName1(e) {
    mergeForm({ signatureName1: e });
  }

  function handleSignatureName2(e) {
    mergeForm({ signatureName2: e });
  }

  return (
    <Form className="p-3">
      <LoadingOverlay active={isFetching} />
      <Row>
        <Col lg="auto">
          <Select
            name="signatureName1"
            label="Signature Name 1"
            value={form.signatureName1}
            disabled={!data}
            options={data}
            onChange={handleSignatureName1}
            control={control}
          />
        </Col>
        <Col lg="auto">
          <Select
            name="signatureName2"
            label="Signature Name 2"
            value={form.signatureName2}
            disabled={!data}
            options={data}
            onChange={handleSignatureName2}
            control={control}
          />
        </Col>
        <Col lg="auto">
          <Group controlId="both">
            <Check
              id="both"
              type="checkbox"
              label="Samples Detected Both Signatures"
              checked={both}
              onChange={(e) => mergeForm({ both: e.target.checked })}
            />
          </Group>
        </Col>
      </Row>
    </Form>
  );
}
