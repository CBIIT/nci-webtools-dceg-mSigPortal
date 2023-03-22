import { Container } from 'react-bootstrap';
import MutProfilesForm from './mutProfiles-form';
import MutProfilesPlot from './mutProfiles-plot';
import { useState } from 'react';

export default function MutationalProfiles({ state }) {
  const [form, setForm] = useState({
    sample: '',
    profile: '',
    matrix: '',
    filter: '',
  });

  const mergeForm = (update) => {
    setForm({ ...form, ...update });
  };

  return (
    <Container
      fluid
      className="bg-white border rounded p-0"
      style={{ minHeight: 500 }}
    >
      <MutProfilesForm state={state} form={form} mergeForm={mergeForm} />
      <hr />
      <MutProfilesPlot state={state} form={form} />
    </Container>
  );
}
