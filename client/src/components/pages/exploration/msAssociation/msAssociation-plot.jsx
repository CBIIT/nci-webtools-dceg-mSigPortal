import React, { useEffect, useState } from 'react';
import {
  useMsAssociationQuery,
  useMsAssociation2SourceQuery,
} from './apiSlice';
import { LoadingOverlay } from '@/components/controls/loading-overlay/loading-overlay';
import Plotly from '@/components/controls/plotly/plot/plot';
import { getDatasetLabel } from '@/services/utils';

export default function MsAssociationPlot({
  state,
  form,
  availableOptions,
  divId = 'msAssociation',
}) {
  const [params, setParams] = useState('');
  const [params2, setParams2] = useState('');
  const { signatureName1, signatureName2, both } = form;
  const { study, strategy, signatureSetName, cancer, useAllCancer, id, id2 } =
    state;

  const { currentData, error, isFetching } = useMsAssociationQuery(params, {
    skip: !params,
  });
  const { data: currentData2, error: error2 } = useMsAssociation2SourceQuery(params2, {
    skip: !params2,
  });

  const missingSignature =
    availableOptions &&
    [signatureName1, signatureName2].find(
      (sig) => sig && !availableOptions.some((e) => e.value === sig.value)
    );

  useEffect(() => {
    if (missingSignature) return;
    if (signatureName1 && signatureName2 && id && !id2) {
      setParams({
        signatureName: signatureName1.value + ';' + signatureName2.value,
        both,
        userId: id,
      });
    } else if (signatureName1 && signatureName2 && id && id2) {
      setParams2([
        {
          signatureName: signatureName1.value,
          userId: signatureName1.id,
          both,
        },
        {
          signatureName: signatureName2.value,
          userId: signatureName2.id,
          both,
        },
      ]);
    } else if (signatureName1 && signatureName2 && study) {
      setParams({
        signatureName: signatureName1.value + ';' + signatureName2.value,
        both,
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        ...(!useAllCancer && { cancer: cancer.value }),
      });
    }
  }, [signatureName1, signatureName2, both, missingSignature]);

  const plotError = error || error2 || currentData?.error || currentData2?.error;
  const plotData = currentData || currentData2;
  const label = getDatasetLabel(state);

  return (
    <div>
      <LoadingOverlay active={isFetching} />
      {label && <p className="text-center text-muted mb-0">{label}</p>}
      {missingSignature ? (
        <p className="p-3 text-danger text-center">
          Signature {missingSignature.value} is missing in this dataset
        </p>
      ) : (
        <>
          {plotError && <p className="p-3 text-danger">{plotError}</p>}
          {plotData && !plotError && (
            <Plotly
              className="w-100"
              divId={divId}
              filename={divId}
              data={plotData.traces}
              layout={plotData.layout}
              config={plotData.config}
            />
          )}
        </>
      )}
    </div>
  );
}
