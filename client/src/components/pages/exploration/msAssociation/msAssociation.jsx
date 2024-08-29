import { Container } from 'react-bootstrap';
import MsAssociationForm from './msAssociation-form';
import MsAssociationPlot from './msAssociation-plot';
import Description from '../../../controls/description/description';
import { useState } from 'react';

export default function MsAssociation({ state }) {
  const [form, setForm] = useState({
    signatureName1: '',
    signatureName2: '',
    both: false,
  });

  const mergeForm = (update) => setForm({ ...form, ...update });

  return (
    <Container fluid className="bg-white border rounded p-0">
      <div className="p-3">
        <b>Mutational Signature Association</b>
        <Description
          less="The scatter plot below illustrates the associations between two selected mutational signatures."
          more="On the x-axis is the number of mutations (log 10) assigned to Signature Name 1, and on the y-axis is the number of mutations (log10) assigned to Signature Name 2. Use the parameter from the top panel, [Samples Detected Both Signatures] to remove the samples that do not carry both signatures before running the association analysis. Use the parameter in the left panel [Cancer Type Only] to perform association analyses on the selected cancer type only."
        />
      </div>
      <hr />
      <MsAssociationForm state={state} form={form} mergeForm={mergeForm} />
      <hr />
      
      <MsAssociationPlot state={state} form={form} />
    </Container>
  );
}
