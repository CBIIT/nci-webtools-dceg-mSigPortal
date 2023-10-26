import { Container } from 'react-bootstrap';
import MsIndividualPlot from './msIndividual-plot';
import MsIndividualForm from './msIndividual-form';
import { useState } from 'react';

export default function MsIndividual({ state }) {
  const [form, setForm] = useState({ sample: '' });
  const mergeForm = (update) => setForm({ ...form, ...update });
  return (
    <Container fluid className="bg-white border rounded p-0">
      <div className="p-3">
        <b>Mutational Signature in Individual Sample</b>
        <div>
          The following plot is used to visualize the signature decomposition in
          individual samples. Select the [Sample Name] and click the
          [Recalculate] button to visualize the signature deconvolution of the
          selected sample.
        </div>
      </div>

      <hr />
      <MsIndividualForm state={state} form={form} mergeForm={mergeForm} />
      <hr />
      <MsIndividualPlot state={state} form={form} />
    </Container>
  );
}
