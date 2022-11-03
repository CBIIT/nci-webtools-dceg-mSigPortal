import { useState, useEffect } from 'react';
import Plotly from '../../../controls/plotly/plot/plot';
import { useSelector } from 'react-redux';
import { useTmbPlotQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

import './plot.scss';
export default function MutProfilePlot() {
  const publicForm = useSelector((state) => state.exposure.publicForm);
  const [params, setParams] = useState('');
  const { data, error, isFetching } = useTmbPlotQuery(params, {
    skip: !params,
  });

  useEffect(() => {
    const { study, strategy, signatureSetName, cancer, useAllCancer } =
      publicForm;
    if (study) {
      setParams({
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        ...(!useAllCancer && { cancer: cancer.value }),
      });
    }
  }, [publicForm]);

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
