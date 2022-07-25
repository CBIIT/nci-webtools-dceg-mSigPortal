import { Container } from 'react-bootstrap';
import MsDecomposition from './msPrevalence-plot';
import Description from '../../../controls/description/description';

export default function MutationalProfiles(props) {
  return (
    <Container fluid className="bg-white border rounded p-0" {...props}>
      <div className="p-3">
        <b>Prevalence of Mutational Signature</b>
        <Description
          className="m-0"
          less="The following plot indicates both mutation and sample level prevalence of signatures from the selected Study."
          more="For prevalence by samples, input the [Minimal Number of Mutations Assigned to Each Signature] to set the smallest number of mutations assigned to each signature required for the detection of the mutational signature in each sample."
        />
      </div>
      <hr />
      <MsDecomposition />
    </Container>
  );
}
