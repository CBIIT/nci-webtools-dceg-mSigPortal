import React, { useEffect, useState } from 'react';
import { useSelector } from 'react-redux';
import axios from 'axios';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import SvgContainer from '../../../controls/svgContainer/svgContainer';
import Description from '../../../controls/description/description';
import { useProfilerSummaryQuery } from './apiSlice';

export default function ProfilerSummary() {
  const store = useSelector((state) => state.visualization);

  const { source, matrixList, projectID } = store.main;
  const { study, cancer, strategy } = store.publicForm;

  const [params, setParams] = useState(null);
  const [plotPath, setPath] = useState('');

  const { data, error, isFetching } = useProfilerSummaryQuery(params, {
    skip: !params,
  });

  useEffect(() => {
    // check if profiler summary already exists, else calculate
    if (projectID) {
      plotExists();
    }
  }, [projectID]);

  // set path after calculation
  useEffect(() => {
    if (data?.output.plotPath) setPath(data.output.plotPath);
  }, [data]);

  async function plotExists() {
    try {
      const path = `${projectID}/results/profilerSummary${
        source == 'public' ? 'Public' : ''
      }/profilerSummary.svg`;
      const response = await axios.head('web/results/' + path);
      setPath(path);
    } catch (error) {
      calculateRequest();
    }
  }

  function calculateRequest() {
    const params =
      source == 'user'
        ? {
            fn: 'profilerSummary',
            args: { matrixList: JSON.stringify(matrixList) },
            projectID,
          }
        : {
            fn: 'profilerSummaryPublic',
            args: {
              study: study.value,
              cancerType: cancer.value,
              experimentalStrategy: strategy.value,
            },
            projectID,
          };

    setParams(params);
  }

  return (
    <div className="bg-white border rounded">
      <LoadingOverlay active={isFetching} />
      <div className="p-3">
        <b>Number of Mutations Per Sample with Regard to Mutational Profile</b>
        <Description
          className="m-0"
          less="This plot illustrates the number of mutations in each tumor sample from [Cancer Type] in the selected [Study]."
          more="On the y-axis is the number of mutations in log base 10 scale, and on the x-axis is the sample index for each sample of the selected cancer type (sorted by number of mutations in ascending order). The different colored lines represent different mutational profiles (SBS= single-base substitution, DBS= doublet-base substitution, ID=indel)."
        />
      </div>

      {plotPath && (
        <>
          <hr />
          <SvgContainer
            title="Number of Mutations Per Sample with Regard to Mutational Profile"
            className="p-3"
            downloadName={plotPath.split('/').slice(-1)[0]}
            plotPath={'web/results/' + plotPath}
            height="600px"
          />
        </>
      )}
      {error && (
        <p className="text-center">
          An error has occured. Please check your inputs and try again.
        </p>
      )}
    </div>
  );
}
