import { useState, useEffect } from 'react';
import Plot from 'react-plotly.js';
import { cloneDeep } from 'lodash';
import { useSelector } from 'react-redux';
import { useTmbPlotQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

export default function MutProfilePlot() {
  const publicForm = useSelector((state) => state.exposure.publicForm);
  const [params, setParams] = useState('');
  const { data, error, isFetching } = useTmbPlotQuery(params, {
    skip: !params,
  });

  useEffect(() => {
    const { study, strategy, signatureSetName } = publicForm;
    if (study) {
      setParams({
        study: study.value,
        strategy: strategy.value,
        signatureSet: signatureSetName.value,
      });
    }
  }, [publicForm]);

  return (
    <div>
      <LoadingOverlay active={isFetching} />
      {data && (
        <Plot
          className="w-100"
          style={{ height: '400px' }}
          data={cloneDeep(data.data)}
          layout={cloneDeep(data.layout)}
          config={cloneDeep(data.config)}
          useResizeHandler
        />
      )}
      {error && <p>An error has occured</p>}
    </div>
  );
}
