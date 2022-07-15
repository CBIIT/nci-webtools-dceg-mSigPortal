import { Container } from 'react-bootstrap';
import TmbSigPlot from './msBurden-plot';
import Description from '../../../controls/description/description';

export default function MutationalProfiles(props) {
  return (
    <Container fluid className="bg-white border rounded p-0" {...props}>
      <div className="p-3">
        <b>Mutational Signature Burden Across Cancer Types</b>
        <Description
          className="m-0"
          less="The bar plot below illustrates mutational signature burden across different cancer types with regard to a specific selected signature."
          more="On the y-axis is the number of mutations per Mb (log10) assigned to selected signatures, and the x-axis denotes the sample numbers. The number in green denotes the number of samples for each cancer type, and the number in blue is the number of samples in that cancer type with the selected signature. "
        />
      </div>
      <hr />
      <TmbSigPlot />
    </Container>
  );
}
