import { useEffect, useState } from 'react';
import Plotly from '../../../controls/plotly/plot/plot';
import { useSelector } from 'react-redux';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Description from '../../../controls/description/description';
import { useProfilerSummaryQuery } from './apiSlice';

export default function ProfilerSummary() {
  const store = useSelector((state) => state.visualization);
  const publicForm = store.publicForm;
  const { source, projectID } = store.main;

  const [params, setParams] = useState();

  const { data, error, isFetching } = useProfilerSummaryQuery(params, {
    skip: !params,
  });

  useEffect(() => {
    if (projectID) {
      setParams({
        study: publicForm.study.value,
        cancer: publicForm.cancer.value,
        strategy: publicForm.strategy.value,
        ...(source == 'user' && { userId: projectID }),
      });
    }
  }, [projectID]);

  return (
    <div className="bg-white border rounded">
      <div className="p-3">
        <b>Number of Mutations Per Sample with Regard to Mutational Profile</b>
        <Description
          className="m-0"
          less="This plot illustrates the number of mutations in each tumor sample from [Cancer Type] in the selected [Study]."
          more="On the y-axis is the number of mutations in log base 10 scale, and on the x-axis is the sample index for each sample of the selected cancer type (sorted by number of mutations in ascending order). The different colored lines represent different mutational profiles (SBS= single-base substitution, DBS= doublet-base substitution, ID=indel)."
        />
      </div>
      <hr />
      <LoadingOverlay active={isFetching} />
      {data && (
        <Plotly
          style={{ height: '500px' }}
          className="w-100"
          data={data.traces}
          layout={data.layout}
          config={data.config}
        />
      )}
      {error && (
        <p className="text-center">
          An error has occured. Please check your inputs and try again.
        </p>
      )}
    </div>
  );
}
