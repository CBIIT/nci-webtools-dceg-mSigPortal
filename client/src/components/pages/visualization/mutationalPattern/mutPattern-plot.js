import { useEffect, useState } from 'react';
import Plotly from '../../../controls/plotly/plot/plot';
import { useSelector } from 'react-redux';
import { useMpeaScatterQuery, useMpeaBarQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

export default function MutPatternPlot() {
  const store = useSelector((state) => state.visualization);
  const publicForm = store.publicForm;
  const { source, projectID, matrixList } = store.main;
  const { proportion, pattern } = store.mutationalPattern;

  const [scatterParams, setScatterParams] = useState('');

  const {
    data: scatterData,
    error: scatterError,
    isFetching: fetchingScatter,
  } = useMpeaScatterQuery(scatterParams, {
    skip: !scatterParams,
  });

  const [barParams, setBarParams] = useState('');
  const {
    data: patternData,
    error: patternError,
    isFetching: fetchingPattern,
  } = useMpeaBarQuery(barParams, {
    skip: !barParams,
  });

  // get data on form change
  useEffect(() => {
    const { study, cancer, strategy } = publicForm;
    if (pattern) {
      setScatterParams({
        profile: 'SBS',
        matrix: '96',
        pattern: pattern,
        ...(source == 'public' && {
          study: study.value,
          cancer: cancer.value,
          strategy: strategy.value,
        }),
        ...(source == 'user' && { userId: projectID }),
      });
    }
  }, [publicForm, proportion, pattern]);

  useEffect(() => {
    if (source == 'public') {
      const { study } = publicForm;
      if (study && proportion) {
        setBarParams({
          study: study.value,
          proportion: proportion,
        });
      }
    } else {
      if (proportion && pattern) {
        setBarParams({
          pattern,
          proportion,
          matrixFile: matrixList.filter(
            (e) => e.profile == 'SBS' && e.matrix == '96'
          )[0].Path,
          userId: projectID,
        });
      }
    }
  }, [publicForm, proportion, matrixList]);

  const divId1 = 'mutationalPatternBarlot';
  const divId2 = 'mutationalPatternScatterlot';
  const config = {
    displayModeBar: true,
    responsive: true,
    displaylogo: false,
    toImageButtonOptions: {
      format: 'svg',
      filename: pattern?.value || 'Mutational Pattern',
      scale: 1,
    },
  };

  return (
    <>
      {patternError && (
        <p className="p-3">An error has occured. Please verify your input.</p>
      )}
      <div id="barchart" style={{ minHeight: '200px' }}>
        <LoadingOverlay active={fetchingPattern} />
        {patternData && patternData?.traces && !patternError && (
          <>
            <Plotly
              className="w-100"
              divId={divId1}
              data={patternData.traces}
              layout={patternData.layout}
              config={config}
            />
            <p className="p-3">
              This plot illustrates the frequency by count of each mutational
              pattern in the given study and cancer type or input dataset. The
              y-axis is the frequency of each mutational pattern across all
              samples, and the x-axis includes each of the mutational patterns
              present in the study and cancer type that meet the criteria for
              the minimal proportion of mutations within each mutational
              pattern.
            </p>
          </>
        )}
        {patternData && !patternData?.traces && !fetchingPattern && (
          <div className="p-3 text-center">
            <b>Frequency of Mutational Pattern</b>
            <p>
              No mutational pattern with proportion of mutations large than{' '}
              {proportion}. Try a lower proportion value.
            </p>
          </div>
        )}
      </div>
      <hr />
      <div id="context" style={{ minHeight: '200px' }}>
        <LoadingOverlay active={fetchingScatter} />
        {scatterData && !scatterError && (
          <>
            <Plotly
              className="w-100"
              divId={divId2}
              data={scatterData.traces}
              layout={scatterData.layout}
              config={config}
            />
            <p className="p-3">
              This plot illustrates the mutational pattern context entered
              compared to other contexts with the same SBS mutation for each
              sample. On the y-axis is the other contexts, and on the x-axis is
              the specific mutational pattern context input for the enrichment
              analysis. For some studies including multiple cancer types (such
              as TCGA PanCancer), different colors will be used.
            </p>
          </>
        )}
      </div>
    </>
  );
}
