import { useRecoilState, useRecoilValue } from 'recoil';
import { formState } from './tmb.state';
import { Form, Row, Col } from 'react-bootstrap';
import Select from 'react-select';

export default function TreeLeafForm() {
  const [form, setForm] = useRecoilState(formState);
  const mergeForm = (state) => setForm({ ...form, ...state });

  function handleChange(value, event) {
    let { name } = event;
    mergeForm({ [name]: value });
  }

  const options = [
    {
      label: 'Single',
      value: {
        study: 'PCAWG',
        strategy: 'WGS',
        signature_set: 'COSMIC_v3_Signatures_GRCh37_SBS96',
        cancer: 'Lung-AdenoCA',
      },
    },
    {
      label: 'Multi',
      value: {
        study: 'PCAWG',
        strategy: 'WGS',
        signature_set: 'COSMIC_v3_Signatures_GRCh37_SBS96',
      },
    },
  ];

  return (
    <Form>
      <Row>
        <Col md="auto">
          <Form.Group controlId="profile" className="mb-3">
            <Form.Label>Cancer Types</Form.Label>
            <Select
              name="option"
              // defaultValue={profileOptions[0]}
              options={options}
              onChange={handleChange}
            />
          </Form.Group>
        </Col>
      </Row>
    </Form>
  );
}
