import React, { useEffect, useState } from 'react';
import { useMsAssociationQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plotly from '../../../controls/plotly/plot/plot';

export default function MsAssociationPlot({ state, form }) {
  const [params, setParams] = useState('');
  const { data, error, isFetching } = useMsAssociationQuery(params, {
    skip: !params,
  });

  const { signatureName1, signatureName2, both } = form;

  const { study, strategy, signatureSetName, cancer, useAllCancer, id } = state;

  useEffect(() => {
    if (signatureName1 && signatureName2 && id) {
      setParams({
        signatureName: signatureName1.value + ';' + signatureName2.value,
        both,
        userId: id,
      });
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
  }, [signatureName1, signatureName2, both, id]);

  return (
    <div>
      {/* <h6 className="p-3 text-center">
        <b>Mutational Signature Association</b>
      </h6> */}
      <LoadingOverlay active={isFetching} />
      {error && <p className="p-3 text-danger">{error}</p>}
      {data && (
        <Plotly
          className="w-100"
          data={data.traces}
          layout={data.layout}
          config={data.config}
        />
      )}
    </div>
  );
}
