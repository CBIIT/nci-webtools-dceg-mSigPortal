import { useEffect, useState } from 'react';
import Plotly from '../../../controls/plotly/plot/plot';
import { useSelector } from 'react-redux';
import { useMutationalPatternScatterQuery, usePatternQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import './plot.scss';
export default function MutPatternPlot() {
  const publicForm = useSelector((state) => state.visualization.publicForm);
  const [scatterParams, setScatterParams] = useState('');

  const {
    data: scatterData,
    error: scatterError,
    isFetching: fetchingScatter,
  } = useMutationalPatternScatterQuery(scatterParams, {
    skip: !scatterParams,
  });

  const [patternParams, setPatternParams] = useState('');
  const {
    data: patternData,
    error: patternError,
    isFetching: fetchingPattern,
  } = usePatternQuery(patternParams, {
    skip: !patternParams,
  });

  const store = useSelector((state) => state.visualization);

  const { proportion, pattern } = store.mutationalPattern;

  // get data on form change
  useEffect(() => {
    const { study, cancer, strategy } = publicForm;
    if (study && pattern) {
      setScatterParams({
        study: study.value,
        cancer: cancer.value,
        strategy: strategy.value,
        profile: 'SBS',
        matrix: '96',
        pattern: pattern,
      });
    }
  }, [publicForm, proportion, pattern]);

  useEffect(() => {
    const { study } = publicForm;
    if (study && proportion) {
      setPatternParams({
        study: study.value,
        proportion: parseFloat(proportion),
      });
    }
  }, [publicForm, proportion]);

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
      <LoadingOverlay active={fetchingPattern} />
      {patternError && (
        <p className="p-3">An error has occured. Please verify your input.</p>
      )}
      <div id="barchart">
        {patternData && !patternError ? (
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
        ) : (
          <div className="text-center my-4">No data available</div>
        )}
      </div>
      <div id="context">
        {scatterData && !scatterError ? (
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
        ) : (
          <div className="text-center my-4">No data available</div>
        )}
      </div>
    </>
  );
}
