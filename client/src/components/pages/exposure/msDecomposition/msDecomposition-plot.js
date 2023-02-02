import { useState, useEffect } from 'react';
import Plotly from '../../../controls/plotly/plot/plot';
import { useMsDecompositionQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

export default function MsDecompositionPlot({ state }) {
  const [params, setParams] = useState('');

  const { data, error, isFetching } = useMsDecompositionQuery(params, {
    skip: !params,
  });

  useEffect(() => {
    const { study, strategy, signatureSetName, cancer } = state;
    if (study) {
      setParams({
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        cancer: cancer.value,
      });
    } else if (state.id) {
      setParams({ userId: state.id });
    }
  }, [state]);

  return (
    <>
      <LoadingOverlay active={isFetching} />
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
