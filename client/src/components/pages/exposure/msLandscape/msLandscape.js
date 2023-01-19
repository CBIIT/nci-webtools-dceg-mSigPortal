import { Container } from 'react-bootstrap';
import MsLandscapePlot from './msLandscape-plot';
import Description from '../../../controls/description/description';

export default function MutationalProfiles(props) {
  function handleVariable(args) {
    console.log(args);
    const reader = new FileReader();
  }

  return (
    <Container fluid className="bg-white border rounded p-0" {...props}>
      <div className="p-3">
        <b>Landscape of Mutational Signature Activity</b>
        <Description
          className="m-0"
          less="The following clustering of mutational signatures allows users to investigate the overall landscape of mutational signatures from the selected cancer type."
          more="The combined plot below includes the following sub-plots from top to bottom: 1) Unsupervised clustering of all samples based on signature contributions; 2) Values from uploaded variable data assigned to each sample (optional); 3) Stacked bar plot that contains the number of mutations assigned to each mutational signature in each sample; 4) Cosine similarity between the original mutational pattern and reconstructed mutational pattern for each sample, which indicates the performance of mutational signature deconvolution; 5) Signature contribution plot, which illustrates the contribution of different signatures within a sample. The colors that correspond to each subplot can be found in the legend."
        />
      </div>
      <hr />
      <MsLandscapePlot handleVariable={handleVariable} />
    </Container>
  );
}
