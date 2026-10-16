import { useState, useEffect } from 'react';
import Plotly from '@/components/controls/plotly/plot/plot';
import { useMsPrevelencePlotQuery } from './apiSlice';
import { LoadingOverlay } from '@/components/controls/loading-overlay/loading-overlay';
import { getDatasetLabel } from '@/services/utils';

export default function MsPrevalancePlot({ form, state, divId = 'msPrevalence' }) {
  const [params, setParams] = useState('');
  const { currentData, error, isFetching } = useMsPrevelencePlotQuery(params, {
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

  const label = getDatasetLabel(state);

  return (
    <div>
      <LoadingOverlay active={isFetching} />
      {label && <p className="text-center text-muted mb-0">{label}</p>}
      {error && <p className="p-3 text-danger">Plot is unavailable</p>}
      {currentData && (
        <div>
          <Plotly
            className="w-100"
            divId={divId}
            filename={divId}
            data={currentData.traces}
            layout={currentData.layout}
            config={currentData.config}
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
