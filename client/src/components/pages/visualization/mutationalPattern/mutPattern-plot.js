import { useEffect, useState } from 'react';
import { useSelector } from 'react-redux';
import SvgContainer from '../../../controls/svgContainer/svgContainer';
import { useMutationalPatternQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

export default function MutProfilePlot() {
  const store = useSelector((state) => state.visualization);

  const { proportion, pattern } = store.mutationalPattern;
  const { projectID, source, matrixList } = store.main;

  const [params, setParams] = useState(null);

  const { data, error, isFetching } = useMutationalPatternQuery(params, {
    skip: !params,
  });

  // get data on form change
  useEffect(() => {
    if (proportion && pattern) {
      const params =
        source == 'user'
          ? {
              fn: 'mutationalPattern',
              args: {
                matrixFile: matrixList.filter(
                  (e) => e.profileType == 'SBS' && e.matrixSize == '96'
                )[0].Path,
                proportion: parseFloat(proportion),
                pattern: pattern,
              },
              projectID,
            }
          : {
              fn: 'mutationalPatternPublic',
              args: {
                study: store.publicForm.study.value,
                cancerType: store.publicForm.cancer.value,
                experimentalStrategy: store.publicForm.strategy.value,
                proportion: parseFloat(proportion),
                pattern: pattern,
              },
              projectID,
            };
      setParams(params);
    }
  }, [proportion, pattern]);

  return (
    <>
      <LoadingOverlay active={isFetching} />
      {error && (
        <p className="p-3">An error has occured. Please verify your input.</p>
      )}

      <div id="barchart">
        {data?.output.barPath && (
          <>
            <SvgContainer
              className="p-3"
              downloadName={data.output.barPath.split('/').slice(-1)[0]}
              plotPath={'results/' + data.output.barPath}
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
        {data?.output.plotPath && !data?.output.barPath && (
          <div className="p-3">
            <p>Frequency of Mutational Pattern</p>
            <p>
              No mutational pattern with proportion of mutations large than{' '}
              {proportion}
            </p>
          </div>
        )}
      </div>
      <div id="context">
        {data?.output.plotPath && (
          <>
            <SvgContainer
              className="p-3"
              downloadName={data.output.plotPath.split('/').slice(-1)[0]}
              plotPath={'results/' + data.output.plotPath}
              txtPath={`results/${data.output.txtPath}`}
              title="Proportion of Mutational Pattern Context Compared to Other Contexts with the same SBS Mutation"
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
