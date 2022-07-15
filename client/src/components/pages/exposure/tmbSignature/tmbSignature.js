import { Container } from 'react-bootstrap';
import TmbSigPlot from './tmbSignature-plot';
import Description from '../../../controls/description/description';

export default function MutationalProfiles(props) {
  return (
    <Container fluid className="bg-white border rounded p-0" {...props}>
      <div className="p-3">
        <b>Tumor Mutational Burden</b>
        <Description
          className="m-0"
          less="The bar plot below illustrates tumor mutational burden (number of mutations per Mb) across different cancer types for the selected study."
          more="The y-axis is the number of mutations per Mb (log10), and the x-axis indicates cancer types. The green number is the number of samples for a given cancer type, and the blue number is the number of samples that have mutation data available for that cancer type."
        />
      </div>
      <hr />
      <TmbSigPlot />
    </Container>
  );
}
