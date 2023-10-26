import { Container } from 'react-bootstrap';
import MsDecompositionPlot from './msDecomposition-plot';
import Description from '../../../controls/description/description';
import { NavHashLink } from 'react-router-hash-link';

export default function MsDecomposition({ state }) {
  return (
    <Container fluid className="bg-white border rounded p-0">
      <div className="p-3">
        <b>Evaluating the Performance of Mutational Signature Decomposition</b>
        <Description
          className="m-0"
          less="The distribution plot below illustrates mutational signature decomposition distribution in a selected cancer type (by selecting [Cancer Type Only] on the left panel) or across different cancer types"
          more={
            <>
              Five different methods are used to measure similarities of the
              original mutational profile matrix and the reconstructed
              mutational profile matrix across all samples: cosine similarity,
              100-L1_Norm_%, 100-L2_Norm_%, KL_Divergence, and Pearson
              Correlation. Click{' '}
              <NavHashLink to="/faq#cosine-similarity">here</NavHashLink> for
              details of these evaluation methods.
            </>
          }
        />
      </div>
      <hr />
      <MsDecompositionPlot state={state} />
    </Container>
  );
}
