import React, { useEffect } from 'react';
import { Form, Row, Col } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import Select from '../../../controls/select/selectHookForm';
import { useSeqmatrixOptionsQuery } from '../../../../services/store/rootApi';

export default function ClusteredForm({ state, form, setForm }) {
  const { id } = state;
  const { sample } = form;

  const { data: options } = useSeqmatrixOptionsQuery(
    { userId: id },
    { skip: !id }
  );

  const { control } = useForm();

  // populate controls
  useEffect(() => {
    if (options && !sample) {
      handleSample(sampleOptions[0]);
    }
  }, [options]);

  const sampleOptions = options
    ? [...new Set(options.map((d) => d.sample))]
        .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
        .map((e) => ({
          label: e,
          value: e,
        }))
    : [];

  function handleSample(sample) {
    setForm({ sample });
  }

  return (
    <Form className="p-3">
      <Row>
        <Col lg="auto">
          <Select
            name="sample"
            label="Sample Name"
            value={sample}
            options={sampleOptions}
            control={control}
            onChange={handleSample}
          />
        </Col>
      </Row>
    </Form>
  );
}
