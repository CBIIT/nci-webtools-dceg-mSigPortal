import { Container } from 'react-bootstrap';
import TmbPlot from './tmb-plot';
import Description from '../../../controls/description/description';

export default function TMB({ state }) {
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
      <TmbPlot state={state} />
    </Container>
  );
}
