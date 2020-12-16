import React, { useEffect } from 'react';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchProfilerSummary,
} from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';
import Debug from '../../controls/debug/debug';

export default function ProfilerSummary({ submitR }) {
  const {
    source,
    study,
    studyOptions,
    cancerType,
    pubExperimentalStrategy,
    pDataOptions,
  } = useSelector((state) => state.visualize);
  const { matrixList, svgList, projectID } = useSelector(
    (state) => state.visualizeResults
  );
  const { plotPath, plotURL, err, debugR, loading } = useSelector(
    (state) => state.profilerSummary
  );

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
        if (source == 'user' && matrixList.length) {
          calculateR('profilerSummary', {
            matrixList: JSON.stringify(matrixList),
          });
        } else if (source == 'public' && svgList.length) {
          calculateR('profilerSummaryPublic', {
            study: study,
            cancerType: cancerType,
            experimentalStrategy: pubExperimentalStrategy,
          });
        }
      }
    };

    if (!loading && !plotPath) {
      checkSummary();
    }
  }, [matrixList, svgList]);

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
          dispatchProfilerSummary({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchProfilerSummary({ err: true, plotURL: '' });
    }
  }

  async function calculateR(fn, args) {
    dispatchProfilerSummary({
      loading: true,
      err: false,
      debugR: '',
    });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();

        dispatchProfilerSummary({
          loading: false,
          debugR: err,
        });
      } else {
        const { debugR, output } = await response.json();
        if (Object.keys(output).length) {
          dispatchProfilerSummary({
            debugR: debugR,
            loading: false,
            plotPath: output.plotPath,
          });
          setRPlot(output.plotPath, 'within');
        } else {
          dispatchProfilerSummary({
            debugR: debugR,
            loading: false,
            err: true,
            plotPath: '',
          });
        }
      }
    } catch (err) {
      dispatchError(err);
      dispatchProfilerSummary({ loading: false });
    }
  }

  return (
    <div>
      <LoadingOverlay active={loading} />
      <Plot
        plotName={plotPath.split('/').slice(-1)[0]}
        plotURL={plotURL}
        maxHeight="600px"
      />
      <Debug msg={debugR} />
    </div>
  );
}
