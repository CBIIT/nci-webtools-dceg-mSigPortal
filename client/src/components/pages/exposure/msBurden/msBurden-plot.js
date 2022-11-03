import { useState, useEffect } from 'react';
import Plotly from '../../../controls/plotly/plot/plot';
import { useSelector } from 'react-redux';
import { useMsBurdenQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

import './plot.scss';
export default function MutProfilePlot() {
  const publicForm = useSelector((state) => state.exposure.publicForm);
  const [params, setParams] = useState('');
  const { data, error, isFetching } = useMsBurdenQuery(params, {
    skip: !params,
  });

  const { signatureName } = useSelector((state) => state.exposure.msBurden);
  const { study, strategy, signatureSetName, cancer, useAllCancer } = publicForm;

  useEffect(() => {
    if (signatureName) {
      setParams({
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        signatureName: signatureName.value,
        ...(!useAllCancer && { cancer: cancer.value }),
      });
    }
  }, [signatureName]);

  return (
    <>
      <LoadingOverlay active={isFetching} />
      <br></br>
      <h5 className="d-flex justify-content-center">
        Mutational Signature Burden Across Cancer Types
      </h5>
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
