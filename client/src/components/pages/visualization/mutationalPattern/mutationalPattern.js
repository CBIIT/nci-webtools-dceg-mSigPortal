import { Container } from 'react-bootstrap';
import { NavHashLink } from 'react-router-hash-link';
import Description from '../../../controls/description/description';
import MutationalPatternPlot from './mutPattern-plot';

export default function MutationalPattern(props) {
  return (
    <Container fluid className="bg-white border rounded p-0" {...props}>
      <div className="p-3">
        <Description
          less={
            <p>
              The aim of the mutational pattern enrichment analysis is to
              determine frequency and enrichment of different types of
              mutational patterns. For more information about mutational pattern
              enrichment, click <NavHashLink to="/faq#mpea">here</NavHashLink>.
            </p>
          }
          more={
            <>
              <p>
                <i>Minimal Proportion:</i> For the “Frequency of Mutational
                Pattern” plot, set the minimal proportion of mutational patterns
                identified in all samples from selected or input study. A
                slightly high proportion, such as 0.5, is suggested.
              </p>
              <p>
                <i>Mutational Pattern:</i> For the enrichment plot of
                “Proportion of Mutational Pattern Context Compared to Other
                Contexts with the same SBS Mutation”, select the mutational
                pattern to identify the enrichment of specific mutation context
                as suggested from “Frequency of Mutational Pattern”. The
                mutational pattern supports common nucleotide symbols. Click
                here for common nucleotide symbol information.
              </p>
            </>
          }
        />
      </div>
      <hr />
      <MutationalPatternPlot />
    </Container>
  );
}
