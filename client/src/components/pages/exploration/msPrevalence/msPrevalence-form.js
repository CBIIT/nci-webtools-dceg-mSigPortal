import { Form, Row, Col, Button } from 'react-bootstrap';
import { useForm, Controller } from 'react-hook-form';

const { Group, Label, Control } = Form;

export default function MsPrevalenceForm({ form, setForm }) {
  const {
    control,
    handleSubmit,
    formState: { errors },
  } = useForm({
    defaultValues: form,
  });

  async function onSubmit(data) {
    setForm(data);
  }

  return (
    <div>
      <hr />
      <Form className="p-3" onSubmit={handleSubmit(onSubmit)}>
        <Row>
          <Col lg="auto">
            <Group
              controlId="prevalenceMutations"
              title="Minimum Number of Mutations Within Each Signature"
            >
              <Label>
                Minimum Number of Mutations Assigned to Each Signature
              </Label>
              <Controller
                name="minimum"
                control={control}
                type="number"
                rules={{ required: true }}
                render={({ field }) => (
                  <Control
                    {...field}
                    placeholder="e.g. 100"
                    isInvalid={errors.minimum}
                  />
                )}
              />

              <Form.Control.Feedback type="invalid">
                Enter a numeric minimum value
              </Form.Control.Feedback>
            </Group>
          </Col>
          <Col lg="auto" className="d-flex">
            <Button className="mt-auto mb-3" variant="primary" type="submit">
              Recalculate
            </Button>
          </Col>
        </Row>
      </Form>
    </div>
  );
}
