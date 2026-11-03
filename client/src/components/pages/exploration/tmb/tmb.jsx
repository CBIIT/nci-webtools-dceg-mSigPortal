import { Container, Row, Col } from 'react-bootstrap';
import TmbPlot from './tmb-plot';
import Description from '@/components/controls/description/description';

export default function TMB({ state, compareState }) {
  const columnLayout =
    compareState && !state.useAllCancer && !compareState.useAllCancer;

  return (
    <Container fluid className="bg-white border rounded p-0">
      <div className="p-3">
        <b>Tumor Mutational Burden</b>
        <Description
          className="m-0"
          less="The bar plot below illustrates tumor mutational burden (number of mutations per Mb) across different cancer types for the selected study."
          more="The y-axis is the number of mutations per Mb (log10), and the x-axis indicates cancer types. The green number is the number of samples for a given cancer type, and the blue number is the number of samples that have mutation data available for that cancer type."
        />
      </div>
      <hr />
      {columnLayout ? (
        <Row className="m-0">
          <Col lg="6" className="p-0">
            <TmbPlot state={state} divId="tmb-primary" />
          </Col>
          <Col lg="6" className="p-0">
            <TmbPlot state={compareState} divId="tmb-compare" />
          </Col>
        </Row>
      ) : (
        <>
          <TmbPlot state={state} divId="tmb-primary" />
          {compareState && (
            <>
              <hr />
              <TmbPlot state={compareState} divId="tmb-compare" />
            </>
          )}
        </>
      )}
    </Container>
  );
}
