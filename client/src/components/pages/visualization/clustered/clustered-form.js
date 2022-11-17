import React, { useEffect } from 'react';
import { Form, Row, Col } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import Select from '../../../controls/select/selectForm';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../../services/store/visualization';

const actions = { ...visualizationActions };

export default function ClusteredForm() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.visualization);

  const mergeState = (state) =>
    dispatch(actions.mergeVisualization({ clustered: state }));

  const { matrixData } = store.main;
  const { sample } = store.clustered;

  const { control } = useForm();

  // populate controls
  useEffect(() => {
    if (matrixData.length && !sample) {
      handleSample(sampleOptions[0]);
    }
  }, [matrixData]);

  const sampleOptions = matrixData.length
    ? [...new Set(matrixData.map((d) => d.sample))]
        .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
        .map((e) => ({
          label: e,
          value: e,
        }))
    : [];

  function handleSample(sample) {
    mergeState({ sample });
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
