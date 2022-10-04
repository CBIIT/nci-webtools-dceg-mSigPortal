import { useState, useEffect } from 'react';
import Plotly from '../../../controls/plotly/plot/plot';
import { Container, Row, Col } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { useMsPrevelencePlotQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

import './plot.scss';
export default function MutProfilePlot() {
  const publicForm = useSelector((state) => state.exposure.publicForm);
  const [params, setParams] = useState('');
  const { data, error, isFetching } = useMsPrevelencePlotQuery(params, {
    skip: !params,
  });

  const { minimum } = useSelector((state) => state.exposure.msPrevalence);

  useEffect(() => {
    const { study, strategy, signatureSetName, cancer } = publicForm;
    if (study) {
      setParams({
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        cancer: cancer.value,
        minimum,
      });
    }
  }, [publicForm, minimum]);

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
