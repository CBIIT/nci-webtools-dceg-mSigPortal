import { Container } from 'react-bootstrap';
import MsAssociationForm from './msAssociation-form';
import MsAssociationPlot from './msAssociation-plot';
import Description from '@/components/controls/description/description';
import { useState, useMemo } from 'react';
import { useMsAssociationOptionsQuery } from './apiSlice';

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

export default function MsAssociation({ state, compareState }) {
  const [form, setForm] = useState({
    signatureName1: '',
    signatureName2: '',
    both: false,
  });

  const mergeForm = (update) => setForm({ ...form, ...update });

  const { study, strategy, signatureSetName, cancer, useAllCancer, id } = state;
  const params = useMemo(
    () => toParams({ study, strategy, signatureSetName, cancer, useAllCancer, id }),
    [study, strategy, signatureSetName, cancer, useAllCancer, id]
  );
  const compareParams = useMemo(
    () => (compareState ? toParams(compareState) : ''),
    [compareState]
  );
  const { data: options } = useMsAssociationOptionsQuery(params, {
    skip: !params,
  });
  const { data: compareOptions } = useMsAssociationOptionsQuery(
    compareParams,
    { skip: !compareParams }
  );

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
      <MsAssociationForm
        state={state}
        compareOptions={compareOptions}
        form={form}
        mergeForm={mergeForm}
      />
      <hr />
      <MsAssociationPlot
        state={state}
        form={form}
        availableOptions={options}
        divId="msAssociation-primary"
      />
      {compareState && (
        <>
          <hr />
          <MsAssociationPlot
            state={compareState}
            form={form}
            availableOptions={compareOptions}
            divId="msAssociation-compare"
          />
        </>
      )}
    </Container>
  );
}
