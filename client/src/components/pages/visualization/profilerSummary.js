import React, { useEffect } from 'react';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';
import Debug from '../../controls/debug/debug';
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
  } = visualization.visualize;
  const { matrixList, projectID } = visualization.results;
  const { filtered } = visualization.mutationalProfiles;
  const {
    plotPath,
    plotURL,
    err,
    debugR,
    loading,
  } = visualization.profilerSummary;
  useEffect(() => {
    // check if profiler summary already exists, else lazy-load calculate
    const checkSummary = async () => {
      const path =
        source == 'user'
          ? `/results/profilerSummary/profilerSummary.svg`
          : `/results/profilerSummaryPublic/profilerSummaryPublic.svg`;
      const check = await fetch(`api/results/${projectID}${path}`, {
        method: 'HEAD',
        cache: 'no-cache',
      });

      if (check.status === 200) {
        setRPlot(path);
      } else {
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
    };

    if (filtered.length && !mainLoading.active && !plotURL) {
      checkSummary();
    }
  }, [filtered, mainLoading]);

  async function setRPlot(plotPath) {
    if (plotPath) {
      try {
        const response = await fetch(`api/results/${projectID}${plotPath}`);
        if (!response.ok) {
          // console.log(await response.json());
        } else {
          const pic = await response.blob();
          const objectURL = URL.createObjectURL(pic);

          if (plotURL) URL.revokeObjectURL(plotURL);
          mergeProfilerSummary({
            plotPath: plotPath,
            plotURL: objectURL,
          });
        }
      } catch (err) {
        mergeError(err.message);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      mergeProfilerSummary({ err: true, plotURL: '' });
    }
  }

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
        const { debugR, output } = await response.json();
        if (Object.keys(output).length) {
          mergeProfilerSummary({
            debugR: debugR,
            loading: false,
          });
          setRPlot(output.plotPath, 'within');
        } else {
          mergeProfilerSummary({
            debugR: debugR,
            loading: false,
            err: true,
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
        <p>
          This plot illustrates the number of mutations in each tumor sample
          from [Cancer Type] in selected [Study]. On the y-axis is the number of
          mutations in log base 10, and on the x-axis is the sample index for
          each sample of the selected cancer type (sorted by number of
          mutations). The legend depicts the different colored lines used on the
          plot to denote different mutational profiles (SBS= single-base
          substitution, DBS= doublet-base substitution, ID=indel).
        </p>
      </div>
      <hr />
      <Plot
        className="p-3"
        downloadNam={plotPath.split('/').slice(-1)[0]}
        plotURL={plotURL}
        maxHeight="600px"
      />
    </div>
  );
}
