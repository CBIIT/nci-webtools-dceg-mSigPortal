import React, { useEffect, useState } from 'react';
import { Form, Row, Col } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import Select from '../../../controls/select/selectForm';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { useMsIndividualOptionQuery } from './apiSlice';


export default function MsIndividualForm({ state, form, mergeForm }) {
  const [signatureOptionParams, setSignatureOptionParams] = useState('');

  const { data, error, isFetching } = useMsIndividualOptionQuery(
    signatureOptionParams,
    {
      skip: !signatureOptionParams,
    }
  );


  const { study, strategy, signatureSetName, cancer, useAllCancer, id } = state;

  const { control } = useForm({
    defaultValues: form,
  });

  // set inital
  useEffect(() => {
    if (data) {
      mergeForm({
        sample: data[0],
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

  function handleSample(e) {
    mergeForm({ sample: e });
  }

  return (
    <Form className="p-3">
      <LoadingOverlay active={isFetching} />
      <Row>
        <Col lg="auto">
          <Select
            name="sampleName"
            label="Sample"
            value={form.sample}
            disabled={!data}
            options={data}
            onChange={handleSample}
            control={control}
          />
        </Col>
      </Row>
    </Form>
  );
}
