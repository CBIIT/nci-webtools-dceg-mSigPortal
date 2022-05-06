import { useRecoilState, useRecoilValue } from 'recoil';
import { formState, getGraphData } from './treeLeaf.state';
import { Form, Row, Col } from 'react-bootstrap';
import MultiSelect from '../../../controls/select/multiSelect';
export default function TreeLeafForm() {
  const { attributes } = useRecoilValue(getGraphData);
  const [form, setForm] = useRecoilState(formState);
  const mergeForm = (state) => setForm({ ...form, ...state });

  function handleChange(event) {
    let { name, value, checked, type } = event.target;
    mergeForm({
      [name]: type === 'checkbox' ? checked : value,
    });
  }

  function handleSearch(e) {
    mergeForm({ searchSamples: e });
  }

  function filterSampleOptions(inputValue = '', limit = 100) {
    return attributes
      .filter(
        (g) =>
          !inputValue ||
          g.Sample.toLowerCase().startsWith(inputValue.toLowerCase())
      )
      .map(({ Sample }) => ({ label: Sample, value: Sample }))
      .slice(0, limit);
  }

  async function handleSearchOptions(inputValue) {
    return filterSampleOptions(inputValue, 40);
  }

  return (
    <Form>
      <Row>
        <Col md="auto">
          <Form.Group controlId="showLabels" className="mb-3">
            <Form.Check
              label="Label Samples"
              type="switch"
              name="showLabels"
              checked={form.showLabels}
              onChange={handleChange}
            />
          </Form.Group>
        </Col>
        <Col md="auto">
          <Form.Group controlId="searchSamples" className="mb-3">
            <Form.Label>Search Samples</Form.Label>
            <MultiSelect
              name="searchSamples"
              placeholder="Sample(s)"
              value={form.searchSamples}
              defaultOptions={filterSampleOptions()}
              loadOptions={handleSearchOptions}
              onChange={handleSearch}
            />
          </Form.Group>
        </Col>
      </Row>
    </Form>
  );
}
