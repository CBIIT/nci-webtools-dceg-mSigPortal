import { useState, useEffect } from 'react';
import Plotly from '../../../controls/plotly/plot/plot';
import { useMsBurdenQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

export default function MsBurdenPlot({ state, form }) {
  const [params, setParams] = useState('');
  const { data, error, isFetching } = useMsBurdenQuery(params, {
    skip: !params,
  });

  const { signatureName } = form;
  const { study, strategy, signatureSetName, cancer, useAllCancer, id } = state;

  // query after signature name is changed
  useEffect(() => {
    if (signatureName && id) {
      setParams({
        signatureName: signatureName.value,
        userId: id,
      });
    } else if (signatureName && study) {
      setParams({
        signatureName: signatureName.value,
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        ...(!useAllCancer && { cancer: cancer.value }),
      });
    }
  }, [signatureName, id]);

  return (
    <>
      <LoadingOverlay active={isFetching} />
      <br></br>
      <h5 className="d-flex justify-content-center">
        Mutational Signature Burden Across Cancer Types
      </h5>
      {data && !error ? (
        <Plotly
          className="w-100"
          data={data.traces}
          layout={data.layout}
          config={data.config}
        />
      ) : (
        <div className="text-center my-4">No data available</div>
      )}
    </>
  );
}
