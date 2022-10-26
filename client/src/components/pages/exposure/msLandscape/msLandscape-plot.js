import { useState, useEffect } from 'react';
import Plotly from '../../../controls/plotly/plot/plot';
import { useSelector } from 'react-redux';
import { useMsLandscapePlotQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

import './plot.scss';
export default function MsLandscapePlot() {
  const publicForm = useSelector((state) => state.exposure.publicForm);
  // const [params, setParams] = useState('');
  // const { data, error, isFetching } = useMsLandscapePlotQuery(params, {
  //   skip: !params,
  // });
  const [calculationQuery, setCalculationQuery] = useState('');
  const { data, error, isFetching } = useMsLandscapePlotQuery(
    calculationQuery,
    {
      skip: !calculationQuery,
    }
  );
  console.log(data);

  // useEffect(() => {
  //   const { study, strategy, signatureSetName } = publicForm;
  //   if (study) {
  //     setParams({
  //       study: study.value,
  //       strategy: strategy.value,
  //       signatureSetName: signatureSetName.value,
  //     });
  //   }
  // }, [publicForm]);
  useEffect(() => {
    const { study, strategy, signatureSetName, cancer } = publicForm;
    if (study) {
      const params_activity = {
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        cancer: cancer.value,
      };
      const params_spectrum = {
        study: study.value,
        strategy: strategy.value,
        cancer: cancer.value,
      };
      const params_signature = {
        signatureSetName: signatureSetName.value,
      };
      setCalculationQuery({
        params_activity,
        params_signature,
        params_spectrum,
      });
    }
  }, [publicForm]);

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
