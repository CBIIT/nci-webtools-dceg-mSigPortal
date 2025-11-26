import { useEffect, useState } from 'react';
import Plotly from '../../../controls/plotly/plot/plot';
import { useMutationalProfilesQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

export default function MutProfilePlot({ state, form }) {
  const { sample, profile, matrix, filter } = form;
  const { study, cancer, strategy, source, id } = state;

  const [params, setParams] = useState(null);

  const { data, error, isFetching } = useMutationalProfilesQuery(params, {
    skip: !params,
  });

  // get data on form change
  useEffect(() => {
    if (sample && sample.value && profile.value && matrix.value) {
      const params = {
        study: study.value,
        cancer: cancer.value,
        strategy: strategy.value,
        ...(source == 'user' && { userId: id }),
        sample: sample.value,
        profile: profile.value,
        matrix: matrix.value,
        ...(filter?.label && { filter: filter.value }),
      };
      setParams(params);
    }
  }, [sample, profile, matrix, filter]);

  return (
    <>
      <LoadingOverlay active={isFetching} />
      {data && !error ? (
        <Plotly
          data={data.traces}
          layout={data.layout}
          config={data.config}
          divId="mutationalProfilePlot"
          filename={sample?.value || 'Mutational Profile'}
        />
      ) : (
        <div className="text-center">
          <div>An error has occurred</div>
          <div>{error?.message}</div>
        </div>
      )}
    </>
  );
}
