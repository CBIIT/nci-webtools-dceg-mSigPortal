import { useState, useEffect } from 'react';
import Plotly from '../../../controls/plotly/plot/plot';
import { useSelector } from 'react-redux';
import { useTmbSignaturesPlotQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

export default function TmbSignaturePlot() {
  const { publicForm, main } = useSelector((state) => state.exposure);
  const [params, setParams] = useState('');
  const { data, error, isFetching } = useTmbSignaturesPlotQuery(params, {
    skip: !params,
  });

  // query after public form is submitted
  useEffect(() => {
    const { study, strategy, signatureSetName, cancer } = publicForm;
    if (study) {
      setParams({
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        cancer: cancer.value,
      });
    }
  }, [publicForm]);
  // query after project id is recieved from user form
  useEffect(() => {
    if (main.id) {
      setParams({ userId: main.id });
    }
  }, [main.id]);

  return (
    <>
      <LoadingOverlay active={isFetching} />
      <br></br>
      <h5 className="d-flex justify-content-center">
        Tumor Mutational Burden Separated by Signatures
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
