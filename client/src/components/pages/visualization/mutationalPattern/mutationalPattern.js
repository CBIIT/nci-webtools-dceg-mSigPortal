import { Container } from 'react-bootstrap';
import MutationalPatternForm from './mutPattern-form';
import MutationalPatternPlot from './mutPattern-plot';

export default function MutationalPattern(props) {
  return (
    <Container fluid className="bg-white border rounded p-0" {...props}>
      <MutationalPatternForm />
      <hr />
      <MutationalPatternPlot />
    </Container>
  );
}
