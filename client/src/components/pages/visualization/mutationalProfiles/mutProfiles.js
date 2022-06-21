import { Container } from 'react-bootstrap';
import MutProfilesForm from './mutProfiles-form';
import MutProfilesPlot from './mutProfiles-plot';

export default function MutationalProfiles(props) {
  return (
    <Container
      fluid
      className="bg-white border rounded p-0"
      style={{ minHeight: 500 }}
      {...props}
    >
      <MutProfilesForm />
      <hr />
      <MutProfilesPlot />
    </Container>
  );
}
