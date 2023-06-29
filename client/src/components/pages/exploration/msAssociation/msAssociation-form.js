import React, { useEffect } from 'react';
import { Form, Row, Col } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import Select from '../../../controls/select/selectHookForm';
import { useMsAssociationOptionsQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

const { Group, Check } = Form;

export default function MsAssociationForm({ state, form, mergeForm }) {
  const { study, strategy, signatureSetName, cancer, useAllCancer, id, id2 } =
    state;

  const {
    data: options,
    error,
    isFetching,
  } = useMsAssociationOptionsQuery(
    study
      ? {
          study: study.value,
          strategy: strategy.value,
          signatureSetName: signatureSetName.value,
          ...(!useAllCancer && { cancer: cancer.value }),
        }
      : { userId: id },
    { skip: !study && !id }
  );
  const {
    data: options2,
    error: error2,
    isFetching: fetchingOptions2,
  } = useMsAssociationOptionsQuery({ userId: id2 }, { skip: !id2 });

  const { both } = form;

  const { control } = useForm({
    defaultValues: form,
  });

  // set initial samples
  useEffect(() => {
    // public or single user data
    if (options && !id2) {
      mergeForm({
        signatureName1: options[0],
        signatureName2: options[1] || options[0],
      });
    }
    if (options && options2) {
      mergeForm({
        signatureName1: options[0],
        signatureName2: options2[0],
      });
    }
    // extraction two sources
  }, [options, options2]);

  function handleSignatureName1(e) {
    mergeForm({ signatureName1: e });
  }

  function handleSignatureName2(e) {
    mergeForm({ signatureName2: e });
  }

  const mergedOptions = [...(options || []), ...(options2 || [])].filter(
    (option, index, self) =>
      self.findIndex((o) => o.value === option.value) === index
  );

  // ...

  return (
    <Form className="p-3">
      <LoadingOverlay active={isFetching} />
      <Row>
        <Col lg="auto">
          <Select
            name="signatureName1"
            label="Signature Name 1"
            value={form.signatureName1}
            disabled={!options}
            options={mergedOptions}
            onChange={handleSignatureName1}
            control={control}
          />
        </Col>
        <Col lg="auto">
          <Select
            name="signatureName2"
            label="Signature Name 2"
            value={form.signatureName2}
            disabled={id2 ? !options2 : !options}
            options={mergedOptions}
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

  return (
    <Form className="p-3">
      <LoadingOverlay active={isFetching} />
      <Row>
        <Col lg="auto">
          <Select
            name="signatureName1"
            label="Signature Name 1"
            value={form.signatureName1}
            disabled={!options}
            //options={[...(options || []), ...(options2 || [])]}
            options={mergedOptions}
            onChange={handleSignatureName1}
            control={control}
          />
        </Col>
        <Col lg="auto">
          <Select
            name="signatureName2"
            label="Signature Name 2"
            value={form.signatureName2}
            disabled={id2 ? !options2 : !options}
            //options={[...(options || []), ...(options2 || [])]}
            options={mergedOptions}
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
