import { useState, useEffect } from 'react';
import Plotly from '../../../controls/plotly/plot/plot';
import { useTmbSignaturesPlotQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

export default function TmbSignaturePlot({ state }) {
  const [params, setParams] = useState('');
  const { currentData, error, isFetching } = useTmbSignaturesPlotQuery(params, {
    skip: !params,
  });
  const { study, strategy, signatureSetName, cancer, id } = state;

  // query after public form is submitted
  useEffect(() => {
    if (study) {
      setParams({
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        cancer: cancer.value,
      });
    }
  }, [study]);
  // query after project id is received from user form
  useEffect(() => {
    if (id) setParams({ userId: id });
  }, [id]);

  return (
    <>
      <LoadingOverlay active={isFetching} />
      <br></br>
      <h5 className="d-flex justify-content-center">
        Tumor Mutational Burden Separated by Signatures
      </h5>
      {currentData && !error ? (
        <Plotly
          className="w-100"
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
