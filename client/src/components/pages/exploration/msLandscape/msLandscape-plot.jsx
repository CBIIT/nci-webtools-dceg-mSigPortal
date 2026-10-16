import { useState, useEffect } from 'react';
import Plotly from '@/components/controls/plotly/plot/plot';
import { useMsLandscapePlotQuery } from './apiSlice';
import { LoadingOverlay } from '@/components/controls/loading-overlay/loading-overlay';
import { readFile, parseMatrix } from '@/components/controls/utils/utils';
import { getDatasetLabel } from '@/services/utils';

export default function MsLandscapePlot({
  state,
  variableFile,
  divId = 'msLandscape',
}) {
  const [params, setParams] = useState('');
  let { currentData, error, isFetching } = useMsLandscapePlotQuery(params, {
    skip: !params,
  });

  async function handleCalculate(state, file) {
    const variableData = file ? parseMatrix(await readFile(file)) : '';
    const { study, strategy, signatureSetName, cancer, id } = state;
    if (study) {
      setParams({
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        cancer: cancer.value,
        variableData,
      });
    } else if (id) {
      setParams({ userId: id, variableData });
    }
  }

  useEffect(() => {
    handleCalculate(state, variableFile);
  }, [state, variableFile]);

  const label = getDatasetLabel(state);

  return (
    <div style={{ minHeight: '500px' }}>
      <LoadingOverlay active={isFetching} />
      {label && <p className="text-center text-muted mb-0">{label}</p>}
      {error && (
        <p className="p-3 text-danger">
          {error.data?.error || error.error || 'Plot is unavailable'}
        </p>
      )}
      {currentData && (
        <Plotly
          divId={divId}
          filename={divId}
          data={currentData.traces}
          layout={currentData.layout}
          config={currentData.config}
        />
      )}
    </div>
  );
}
