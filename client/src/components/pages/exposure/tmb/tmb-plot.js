import { useState, useEffect } from 'react';
import Plotly from '../../../controls/plotly/plot/plot';
import { useTmbPlotQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

export default function TmbPlot({ state }) {
  const [params, setParams] = useState('');
  const { data, error, isFetching } = useTmbPlotQuery(params, {
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

  return (
    <>
      <LoadingOverlay active={isFetching} />
      <br></br>
      <h5 className="d-flex justify-content-center">Tumor Mutational Burden</h5>
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
