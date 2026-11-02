import { useState, useEffect } from 'react';
import Plotly from '@/components/controls/plotly/plot/plot';
import { useMsDecompositionQuery } from './apiSlice';
import { LoadingOverlay } from '@/components/controls/loading-overlay/loading-overlay';
import { getDatasetLabel } from '@/services/utils';

export default function MsDecompositionPlot({ state, divId = 'msDecomposition' }) {
  const [params, setParams] = useState('');

  const { currentData, error, isFetching } = useMsDecompositionQuery(params, {
    skip: !params,
  });

  const { study, strategy, signatureSetName, cancer, id } = state;
  useEffect(() => {
    if (study) {
      setParams({
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        cancer: cancer.value,
      });
    } else if (id) {
      setParams({ userId: id });
    }
  }, [study, id]);

  const label = getDatasetLabel(state);

  return (
    <>
      <LoadingOverlay active={isFetching} />
      {label && <p className="text-center text-muted mb-0">{label}</p>}
      {error && <p className="p-3 text-danger">Plot is unavailable</p>}
      {currentData && !error ? (
        <Plotly
          className="w-100"
          divId={divId}
          filename={divId}
          data={currentData.traces}
          layout={currentData.layout}
          config={currentData.config}
        />
      ) : (
        <div className="text-center my-4">No data available</div>
      )}
    </>
  );
}
