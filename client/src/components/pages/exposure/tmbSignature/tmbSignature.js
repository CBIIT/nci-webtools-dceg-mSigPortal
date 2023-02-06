import { Container } from 'react-bootstrap';
import TmbSigPlot from './tmbSignature-plot';
import Description from '../../../controls/description/description';

export default function TmbSignature({ state }) {
  return (
    <Container fluid className="bg-white border rounded p-0">
      <div className="p-3">
        <b>Tumor Mutational Burden Separated by Signatures</b>
        <Description
          className="m-0"
          less="The bar plot below illustrates tumor mutational burden across different signatures from the selected study and cancer type."
          more="On the y-axis is the number of mutations per Mb (log10), and the x-axis denotes sample numbers. The green number is the total number of samples with the selected cancer type (and therefore evaluated for each signature). The blue number is the number of samples with the selected cancer type that detected the mutational signature from the reference signature set."
        />
      </div>
      <hr />
      <TmbSigPlot state={state} />
    </Container>
  );
}
