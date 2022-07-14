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
  console.log(data);
  console.log(publicForm);
  if (data) {
    console.log(data.traces.length);
  }

  return (
    <div
      class="container d-flex align-items-center justify-content-center"
      style={{ minHeight: '500px' }}
    >
      <LoadingOverlay active={isFetching} />
      {data && (
        <Plot
          {...(data.traces.length > 1
            ? { className: 'w-100' }
            : { className: 'w-30' })}
          style={{ height: '500px' }}
          data={cloneDeep(data.traces)}
          layout={cloneDeep(data.layout)}
          config={cloneDeep(data.config)}
          useResizeHandler
        />
      )}
      {error && <p>An error has occured</p>}
    </div>
  );
}
