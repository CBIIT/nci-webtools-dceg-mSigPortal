import React, { useEffect } from 'react';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../services/store/visualization';
import { actions as modalActions } from '../../../services/store/modal';
import Description from '../../controls/description/description';

const actions = { ...visualizationActions, ...modalActions };
export default function ProfilerSummary({ submitR }) {
  const dispatch = useDispatch();
  const visualization = useSelector((state) => state.visualization);
  const mergeProfilerSummary = (state) =>
    dispatch(actions.mergeVisualization({ profilerSummary: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));
  const {
    source,
    study,
    cancerType,
    pubExperimentalStrategy,
    loading: mainLoading,
    matrixList,
    projectID,
  } = visualization.state;

  const { plotPath, err, debugR, loading } = visualization.profilerSummary;
  useEffect(() => {
    // check if profiler summary already exists, else lazy-load calculate
    if (!plotPath && !mainLoading.active && !loading && projectID) {
      if (source == 'user') {
        calculateR('profilerSummary', {
          matrixList: JSON.stringify(matrixList),
        });
      } else if (source == 'public') {
        calculateR('profilerSummaryPublic', {
          study: study,
          cancerType: cancerType,
          experimentalStrategy: pubExperimentalStrategy,
        });
      }
    }
  }, [mainLoading.active, projectID]);

  async function calculateR(fn, args) {
    mergeProfilerSummary({
      loading: true,
      err: false,
      debugR: '',
    });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();

        mergeProfilerSummary({
          loading: false,
          debugR: err,
        });
      } else {
        const { output } = await response.json();
        if (output.plotPath) {
          mergeProfilerSummary({
            debugR: debugR,
            loading: false,
            plotPath: output.plotPath,
          });
        } else {
          mergeProfilerSummary({
            loading: false,
            err: output.error || output.uncaughtError || true,
            plotPath: '',
          });
        }
      }
    } catch (err) {
      mergeError(err.message);
      mergeProfilerSummary({ loading: false });
    }
  }

  return (
    <div className="bg-white border rounded">
      <LoadingOverlay active={loading} />
      <div className="p-3">
        <b>Number of Mutations Per Sample with Regard to Mutational Profile</b>
        <Description
          className="m-0"
          less="This plot illustrates the number of mutations in each tumor sample from [Cancer Type] in the selected [Study]."
          more="On the y-axis is the number of mutations in log base 10 scale, and on the x-axis is the sample index for each sample of the selected cancer type (sorted by number of mutations in ascending order). The different colored lines represent different mutational profiles (SBS= single-base substitution, DBS= doublet-base substitution, ID=indel)."
        />
        {/* <button
          onClick={() =>
            source == 'user'
              ? calculateR('profilerSummary', {
                  matrixList: JSON.stringify(matrixList),
                })
              : calculateR('profilerSummaryPublic', {
                  study: study,
                  cancerType: cancerType,
                  experimentalStrategy: pubExperimentalStrategy,
                })
          }
        >
          calculate
        </button> */}
      </div>
      <div>
        {plotPath && (
          <>
            <hr />
            <Plot
              title="Number of Mutations Per Sample with Regard to Mutational Profile"
              className="p-3"
              downloadName={plotPath.split('/').slice(-1)[0]}
              plotPath={'web/results/' + plotPath}
              height="600px"
            />
          </>
        )}
      </div>
    </div>
  );
}
