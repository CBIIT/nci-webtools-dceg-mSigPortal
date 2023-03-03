import { useState, useEffect } from 'react';
import Plotly from '../../../controls/plotly/plot/plot';

import { useMsLandscapePlotQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { readFile, parseMatrix } from '../../../controls/utils/utils';

export default function MsLandscapePlot({ state, variableFile }) {
  const [params, setParams] = useState('');
  let { data, error, isFetching } = useMsLandscapePlotQuery(params, {
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

  return (
    <div style={{ minHeight: '500px' }}>
      <LoadingOverlay active={isFetching} />
      {error && <p className="p-3 text-danger">{error}</p>}
      {data && (
        <Plotly data={data.traces} layout={data.layout} config={data.config} />
      )}
    </div>
  );
}
