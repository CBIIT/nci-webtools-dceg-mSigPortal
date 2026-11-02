import React, { useEffect } from 'react';
import { Form, Row, Col } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import Select from '@/components/controls/select/selectHookForm';

export default function MsBurdenForm({ options, compareOptions, form, setForm }) {
  const { control } = useForm({ defaultValues: form });

  const mergedOptions = [...(options || []), ...(compareOptions || [])]
    .filter(
      (option, index, self) =>
        self.findIndex((o) => o.value === option.value) === index
    )
    .sort((a, b) =>
      a.label.localeCompare(b.label, undefined, { numeric: true })
    );

  // set inital
  useEffect(() => {
    if (mergedOptions.length) {
      setForm({ signatureName: mergedOptions[0] });
    }
  }, [options, compareOptions]);

  return (
    <Form className="p-3">
      <Row>
        <Col lg="auto">
          <Select
            name="signatureName"
            label="Signature Name"
            value={form.signatureName}
            disabled={!options}
            control={control}
            options={mergedOptions}
            onChange={(name) => setForm({ signatureName: name })}
          />
        </Col>
      </Row>
    </Form>
  );
}
