import { useEffect, useState } from 'react';
import { Container, Row, Col } from 'react-bootstrap';
import Plot from 'react-plotly.js';
import { useSelector } from 'react-redux';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Description from '../../../controls/description/description';
import profilerSummary from '../../../controls/plotly/profilerSummary/profilerSummary';

export default function ProfilerSummary() {
  const store = useSelector((state) => state.visualization);
  const { matrixData } = store.main;

  const [plot, setPlot] = useState({});

  useEffect(() => {
    if (matrixData.length) {
      const plot = profilerSummary(matrixData);
      setPlot(plot);
    }
  }, [matrixData]);

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

      <Container fluid style={{ minHeight: '500px' }} className="mb-3">
        <Row>
          <Col>
            {plot ? (
              <Plot
                className="w-100"
                data={plot.traces}
                layout={plot.layout}
                config={plot.config}
                useResizeHandler
              />
            ) : (
              <LoadingOverlay active={true} />
            )}
          </Col>
        </Row>
      </Container>
      {false && (
        <p className="text-center">
          An error has occured. Please check your inputs and try again.
        </p>
      )}
    </div>
  );
}
