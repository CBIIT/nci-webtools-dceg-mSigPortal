import { Container, Row, Col } from 'react-bootstrap';
import TmbBurdenPlot from './msBurden-plot';
import TmbBurdenForm from './msBurden-form';
import Description from '@/components/controls/description/description';
import { useState, useEffect } from 'react';
import { useMsBurdenOptionsQuery } from './apiSlice';

function toParams(state) {
  if (!state) return '';
  const { study, strategy, signatureSetName, cancer, useAllCancer, id } = state;
  if (study) {
    return {
      study: study.value,
      strategy: strategy.value,
      signatureSetName: signatureSetName.value,
      ...(!useAllCancer && { cancer: cancer.value }),
    };
  } else if (id) {
    return { userId: id };
  }
  return '';
}

export default function MsBurden({ state, compareState }) {
  const [form, setForm] = useState({ signatureName: '' });
  const [params, setParams] = useState('');
  const [compareParams, setCompareParams] = useState('');

  const { data: options } = useMsBurdenOptionsQuery(params, { skip: !params });
  const { data: compareOptions } = useMsBurdenOptionsQuery(compareParams, {
    skip: !compareParams,
  });

  useEffect(() => setParams(toParams(state)), [state]);
  useEffect(
    () => setCompareParams(compareState ? toParams(compareState) : ''),
    [compareState]
  );

  const columnLayout =
    compareState && !state.useAllCancer && !compareState.useAllCancer;

  const primaryPlot = (
    <TmbBurdenPlot
      form={form}
      state={state}
      availableOptions={options}
      divId="msBurden-primary"
    />
  );
  const comparePlot = compareState && (
    <TmbBurdenPlot
      form={form}
      state={compareState}
      availableOptions={compareOptions}
      divId="msBurden-compare"
    />
  );

  return (
    <Container fluid className="bg-white border rounded p-0">
      <div className="p-3">
        <b>Mutational Signature Burden Across Cancer Types</b>
        <Description
          className="m-0"
          less="The bar plot below illustrates mutational signature burden across different cancer types with regard to a specific selected signature."
          more="On the y-axis is the number of mutations per Mb (log10) assigned to selected signatures, and the x-axis denotes the sample numbers. The number in green denotes the number of samples for each cancer type, and the number in blue is the number of samples in that cancer type with the selected signature. "
        />
      </div>
      <hr />
      <TmbBurdenForm
        options={options}
        compareOptions={compareOptions}
        form={form}
        setForm={setForm}
      />
      <hr />
      {columnLayout ? (
        <Row className="m-0">
          <Col lg="6" className="p-0">
            {primaryPlot}
          </Col>
          <Col lg="6" className="p-0">
            {comparePlot}
          </Col>
        </Row>
      ) : (
        <>
          {primaryPlot}
          {comparePlot && (
            <>
              <hr />
              {comparePlot}
            </>
          )}
        </>
      )}
    </Container>
  );
}
