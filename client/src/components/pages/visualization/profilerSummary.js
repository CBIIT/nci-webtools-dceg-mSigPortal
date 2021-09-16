import React, { useEffect } from 'react';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../services/store/visualization';
import { actions as modalActions } from '../../../services/store/modal';

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
      <div className="p-3">
        <b>Number of Mutations Per Sample with Regard to Mutational Profile</b>
        <p>
          This plot illustrates the number of mutations in each tumor sample
          {source == 'public' &&
            ` from Cancer Type: ${cancerType} in selected Study: ${study}`}
          . On the y-axis is the number of mutations in log base 10, and on the
          x-axis is the sample index for each sample of the selected cancer type
          (sorted by number of mutations). The legend depicts the different
          colored lines used on the plot to denote different mutational profiles
          (SBS= single-base substitution, DBS= doublet-base substitution,
          ID=indel).
        </p>
      </div>
      <hr />
      <div style={{ minHeight: '500px' }}>
        <LoadingOverlay active={loading} />
        {plotPath && (
          <Plot
            className="p-3"
            downloadName={plotPath.split('/').slice(-1)[0]}
            plotPath={'api/results/' + plotPath}
            height="600px"
          />
        )}
      </div>
    </div>
  );
}
