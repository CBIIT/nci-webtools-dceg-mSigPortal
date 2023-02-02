import { useState, useEffect } from 'react';
import Plotly from '../../../controls/plotly/plot/plot';
import { useSelector } from 'react-redux';
import { useMsDecompositionQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

export default function MsDecompositionPlot() {
  const { publicForm, main } = useSelector((state) => state.exposure);
  const [params, setParams] = useState('');
  // const { data, error, isFetching } = useMsDecompositionQuery(params, {
  //   skip: !params,
  // });
  // const [calculationQuery, setCalculationQuery] = useState('');
  const { data, error, isFetching } = useMsDecompositionQuery(params, {
    skip: !params,
  });

  console.log(data);
  useEffect(() => {
    const { study, strategy, signatureSetName, cancer } = publicForm;
    if (study) {
      setParams({
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        cancer: cancer.value,
      });
    } else if (main.id) {
      setParams({ userId: main.id });
    }
  }, [publicForm, main]);

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
