import { useEffect, useState } from 'react';
import cloneDeep from 'lodash/cloneDeep';
import { Button, Container, Row, Col } from 'react-bootstrap';
import Plot from 'react-plotly.js';
import { useSelector } from 'react-redux';
import SvgContainer from '../../../controls/svgContainer/svgContainer';
import { useMutationalPattern2Query } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import './plot.scss';
export default function MutProfilePlot() {
  const publicForm = useSelector((state) => state.exposure.publicForm);
  const [params, setParams] = useState('');
  const { data, error, isFetching } = useMutationalPattern2Query(params, {
    skip: !params,
  });
  const store = useSelector((state) => state.visualization);

  const { proportion, pattern } = store.mutationalPattern;
  const { projectID, source, matrixList } = store.main;

  // get data on form change

  useEffect(() => {
    const { study, cancer, strategy, proportion, pattern } = publicForm;
    if (study) {
      setParams({
        study: study.value,
        cancerType: cancer.value,
        experimentalStrategy: strategy.value,
        proportion: parseFloat(proportion),
        pattern: pattern,
      });
    }
  }, [proportion, pattern]);

  const divId = 'mutationalPatternlot';
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
      <LoadingOverlay active={isFetching} />
      {error && (
        <p className="p-3">An error has occured. Please verify your input.</p>
      )}

      <div id="barchart">
        {data && (
          <>
            <Container fluid style={{ minHeight: '500px' }} className="mb-3">
              <Row>
                <Col>
                  <Plot
                    className="w-100"
                    divId={divId}
                    data={cloneDeep(data.traces)}
                    layout={cloneDeep(data.layout)}
                    config={cloneDeep(config)}
                    useResizeHandler
                  />
                </Col>
              </Row>
            </Container>
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
              plotPath={'web/results/' + data.output.plotPath}
              txtPath={`web/results/${data.output.txtPath}`}
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
