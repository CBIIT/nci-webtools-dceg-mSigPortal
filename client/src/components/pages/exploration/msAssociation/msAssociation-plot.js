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

  const { data, error, isFetching } = useMsAssociationQuery(params, {
    skip: !params,
  });
  const {
    data: alt,
    error: altError,
    isFetching: altFetching,
  } = useMsAssociation2SourceQuery(params2, {
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

  return (
    <div>
      <LoadingOverlay active={isFetching} />
      {(error || altError) && (
        <p className="p-3 text-danger">{error || altError}</p>
      )}
      {(data || alt) && (
        <Plotly
          className="w-100"
          data={data?.traces || alt?.traces}
          layout={data?.layout || alt?.layout}
          config={data?.config || alt?.config}
        />
      )}
    </div>
  );
}
