import { useEffect, useState } from 'react';

import Plotly from '../../../controls/plotly/plot/plot';
import { useSelector } from 'react-redux';
import { useMutationalProfilesQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

export default function MutProfilePlot() {
  const store = useSelector((state) => state.visualization);

  const { sample, profile, matrix } = store.mutationalProfiles;
  const { study, cancer, strategy } = store.publicForm;
  const { source, projectID } = store.main;

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
        ...(source == 'user' && { userId: projectID }),
        sample: sample.value,
        profile: profile.value,
        matrix: matrix.value,
      };
      setParams(params);
    }
  }, [sample, profile, matrix]);

  return (
    <>
      <LoadingOverlay active={isFetching} />
      {error ? (
        <div className="text-center">
          <div>An error as occured</div>
          <div>{error.message}</div>
        </div>
      ) : (
        data && (
          <Plotly
            data={data.traces}
            layout={data.layout}
            config={data.config}
            divId="mutationalProfilePlot"
            filename={sample?.value || 'Mutational Profile'}
          />
        )
      )}
    </>
  );
}
