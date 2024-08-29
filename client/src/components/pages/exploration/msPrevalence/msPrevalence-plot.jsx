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
    <div>
      <LoadingOverlay active={isFetching} />
      {error && <p className="p-3 text-danger">Plot is unavailable</p>}
      {data && (
        <div>
          <Plotly
            className="w-100"
            data={data.traces}
            layout={data.layout}
            config={data.config}
          />
          <div className="p-3">
            <p>
              The pie chart on the left illustrates the prevalence of each
              mutational signature by mutations. The bar plot on the right
              illustrates the prevalence of each mutational signature by
              samples. The colors represent the mutational signatures in both
              plots.
            </p>
          </div>
        </div>
      )}
    </div>
  );
}
