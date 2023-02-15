import React, { useEffect, useState } from 'react';
import { Form, Row, Col } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import Select from '../../../controls/select/selectForm';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { useMsIndividualOptionQuery } from './apiSlice';

const { Group, Check } = Form;

export default function MsIndividualForm({ state, form, mergeForm }) {
  const [signatureOptionParams, setSignatureOptionParams] = useState('');

  const { data, error, isFetching } = useMsIndividualOptionQuery(
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
        sampleName: data[0],
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

  function handleSampleName(e) {
    mergeForm({ sampleName: e });
  }

  return (
    <Form className="p-3">
      <LoadingOverlay active={isFetching} />
      <Row>
        <Col lg="auto">
          <Select
            name="sampleName"
            label="Sample"
            value={form.sampleName}
            disabled={!data}
            options={data}
            onChange={handleSampleName}
            control={control}
          />
        </Col>
      </Row>
    </Form>
  );
}
