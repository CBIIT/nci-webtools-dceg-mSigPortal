import { Container } from 'react-bootstrap';
import Description from '../../../controls/description/description';
import { NavHashLink } from 'react-router-hash-link';
import MsIndividualPlot from './msIndividual-plot';
import MsIndividualForm from './msIndividual-form';
import { useState } from 'react';

export default function MsIndividual({ state }) {
  const [form, setForm] = useState({ sample: '' });
  const mergeForm = (update) => setForm({ ...form, ...update });
  return (
    <Container fluid className="bg-white border rounded p-0">
      <div className="p-3">
        <b>Mutational Signature in Individual Sample</b>
        <Description
          less="The following plot is used to visualize the signature decomposition in individual samples. Select the [Sample Name] and click the [Recalculate] button to
          visualize the signature deconvolution of the selected sample."
          more={
            <>
              <p className="mt-3">
                The combination plot shows the original mutational profile, the
                deconstructed mutational profile, and the difference of each
                mutation type between these two profiles, mutational signature
                profiles and proportion of each contributed signature detected
                in the selected sample. Two measurements (RSS and Cosine
                Similarity) of evaluating signature deconvolution are shown on
                the top of this plot.¸
              </p>
              <p>
                Residual Sum of Squares (RSS) measures the discrepancy between
                two mutational profiles. Cosine similarity measures how similar
                two mutational profiles are. For example, two identical
                mutational profiles will have RSS = 0 and Cosine similarity = 1.
                For additional information about RSS and cosine similarity,
                click{' '}
                <NavHashLink to="/faq#cosine-similarity">here</NavHashLink>.
              </p>
            </>
          }
        />
      </div>
      <hr />
      <MsIndividualForm state={state} form={form} mergeForm={mergeForm} />
      <hr />
      <MsIndividualPlot state={state} form={form} setForm={setForm} />
    </Container>
  );
}
