import { useState } from 'react';
import { Container } from 'react-bootstrap';
import Description from '../../../controls/description/description';
import MsPrevalencePlot from './msPrevalence-plot';
import MsPrevalenceForm from './msPrevalence-form';

export default function MsPrevalence({ state }) {
  const [form, setForm] = useState({ minimum: 100 });

  return (
    <Container fluid className="bg-white border rounded p-0">
      <div className="p-3">
        <b>Prevalence of Mutational Signature</b>
        <Description
          className="m-0"
          less="The following plot indicates both mutation and sample level prevalence of signatures from the selected Study."
          more="For prevalence by samples, input the [Minimal Number of Mutations Assigned to Each Signature] to set the smallest number of mutations assigned to each signature required for the detection of the mutational signature in each sample."
        />
      </div>
      <MsPrevalenceForm form={form} setForm={setForm} />
      <hr />

      <MsPrevalencePlot form={form} state={state} />
    </Container>
  );
}
