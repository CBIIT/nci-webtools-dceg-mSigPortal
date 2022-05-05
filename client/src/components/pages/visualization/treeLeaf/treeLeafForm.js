import { useRecoilState } from 'recoil';
import { formState } from './treeLeaf.state';
import { Form, Row, Col } from 'react-bootstrap';

export default function TreeLeafForm() {
  const [form, setForm] = useRecoilState(formState);
  const mergeForm = (state) => setForm({ ...form, ...state });

  function handleChange(event) {
    let { name, value, checked, type } = event.target;
    mergeForm({
      [name]: type === 'checkbox' ? checked : value,
    });
  }

  return (
    <Form>
      <Row>
        <Col>
          <Form.Group controlId="showLabels" className="mb-3">
            <Form.Check
              label="Show Labels"
              type="switch"
              name="showLabels"
              checked={form.showLabels}
              onChange={handleChange}
            />
          </Form.Group>
        </Col>
      </Row>
    </Form>
  );
}
