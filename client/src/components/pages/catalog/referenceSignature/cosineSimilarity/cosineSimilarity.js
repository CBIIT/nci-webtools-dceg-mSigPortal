import { NavHashLink } from 'react-router-hash-link';
import Description from '../../../../controls/description/description';
import CosineSimilarityPlot from './form-plot';

export default function CosineSimilarity() {
  return (
    <div>
      <Description
        className="p-3"
        less="Cosine similarity is a measure of the similarity of two signature matrices, which can be helpful to compare two mutational profiles or signatures. Below you can explore cosine similarity between two reference mutational signature sets."
        more={
          <span>
            Use the dropdown menus to enter a [Profile Name], [Reference
            Signature Set 1], and [Reference Signature Set 2]. Click{' '}
            <NavHashLink to="/faq#cosine-similarity">here</NavHashLink> to learn
            more about cosine similarity.
          </span>
        }
      />
      <hr />
      <CosineSimilarityPlot />
    </div>
  );
}
