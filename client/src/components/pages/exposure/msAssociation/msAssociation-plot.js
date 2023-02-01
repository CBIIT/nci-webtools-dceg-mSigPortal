import { useState, useEffect } from 'react';
import Plotly from '../../../controls/plotly/plot/plot';
import { useSelector } from 'react-redux';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { useMsAssociationOptionsQuery } from './apiSlice';

export default function MsAssociationPlot() {
  const { publicForm, main } = useSelector((state) => state.exposure);
  const [params, setParams] = useState('');
  const { data, error, isFetching } = useMsAssociationOptionsQuery(params, {
    skip: !params,
  });
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
