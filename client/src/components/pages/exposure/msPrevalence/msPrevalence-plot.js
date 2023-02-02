import { useState, useEffect } from 'react';
import Plotly from '../../../controls/plotly/plot/plot';
import { useMsPrevelencePlotQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

export default function MsPrevalancePlot({ form, state }) {
  const [params, setParams] = useState('');
  const { data, error, isFetching } = useMsPrevelencePlotQuery(params, {
    skip: !params,
  });

  const { minimum } = form;

  useEffect(() => {
    const { study, strategy, signatureSetName, cancer } = state;
    if (study) {
      setParams({
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        cancer: cancer.value,
        minimum,
      });
    } else if (state.id) {
      setParams({ userId: state.id, minimum });
    }
  }, [form, state]);

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
