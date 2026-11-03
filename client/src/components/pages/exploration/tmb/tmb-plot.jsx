import { useState, useEffect } from 'react';
import Plotly from '@/components/controls/plotly/plot/plot';
import { useTmbPlotQuery } from './apiSlice';
import { LoadingOverlay } from '@/components/controls/loading-overlay/loading-overlay';
import { getDatasetLabel } from '@/services/utils';

export default function TmbPlot({ state, divId = 'tmb' }) {
  const [params, setParams] = useState('');
  const { currentData, error, isFetching } = useTmbPlotQuery(params, {
    skip: !params,
  });

  const { study, strategy, signatureSetName, cancer, useAllCancer, id } = state;
  // query after public form is submitted
  useEffect(() => {
    if (study) {
      setParams({
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        ...(!useAllCancer && { cancer: cancer.value }),
      });
    }
  }, [study]);
  // query after project id is recieved from user form
  useEffect(() => {
    if (id) {
      setParams({ userId: id });
    }
  }, [id]);

  const label = getDatasetLabel(state);

  return (
    <>
      <LoadingOverlay active={isFetching} />
      <br></br>
      <h5 className="d-flex justify-content-center">Tumor Mutational Burden</h5>
      {label && (
        <p className="text-center text-muted mb-0">{label}</p>
      )}
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
