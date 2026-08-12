import React, { useEffect, useState } from 'react';
import {
  useMsAssociationQuery,
  useMsAssociation2SourceQuery,
} from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plotly from '../../../controls/plotly/plot/plot';

export default function MsAssociationPlot({ state, form }) {
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

  useEffect(() => {
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
  }, [signatureName1, signatureName2, both]);

  const plotError = error || error2 || currentData?.error || currentData2?.error;
  const plotData = currentData || currentData2;
  console.log('plotError', plotError);
  console.log('plotData', plotData);
  return (
    <div>
      <LoadingOverlay active={isFetching} />
      {plotError && <p className="p-3 text-danger">{plotError}</p>}
      {plotData && !plotError && (
        <Plotly
          className="w-100"
          data={plotData.traces}
          layout={plotData.layout}
          config={plotData.config}
        />
      )}
    </div>
  );
}
