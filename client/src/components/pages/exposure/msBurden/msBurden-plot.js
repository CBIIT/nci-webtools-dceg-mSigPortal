import { useState, useEffect } from 'react';
import Plotly from '../../../controls/plotly/plot/plot';
import { useSelector } from 'react-redux';
import { useMsBurdenQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

export default function MsBurdenPlot() {
  const { publicForm, main } = useSelector((state) => state.exposure);
  const [params, setParams] = useState('');
  const { data, error, isFetching } = useMsBurdenQuery(params, {
    skip: !params,
  });

  const { signatureName } = useSelector((state) => state.exposure.msBurden);
  const { study, strategy, signatureSetName, cancer, useAllCancer } =
    publicForm;

  // query after signature name is changed
  useEffect(() => {
    if (signatureName && main.id) {
      setParams({
        signatureName: signatureName.value,
        userId: main.id,
      });
    } else if (signatureName && study) {
      setParams({
        signatureName: signatureName.value,
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        ...(!useAllCancer && { cancer: cancer.value }),
      });
    }
  }, [signatureName, main]);

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
