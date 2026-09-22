import { useSelector } from 'react-redux';
import Select from 'react-select';
import { Form, Row, Col } from 'react-bootstrap';
import MultiSelect from '../../../controls/select/multiSelect';
import { colorOptions, userColorOptions } from './treeLeaf.state';

export default function TreeLeafForm({ isUser, form, onChange, attributes }) {
  const publicForm = useSelector((store) => store.visualization.publicForm);
  const cancers = publicForm?.cancers?.filter((c) => c.value !== '*ALL') || [];
  // only show "All" option if the study doesn't already have an "ALL" cancer type
  const cancerTypes = cancers.some((c) => c.label === 'ALL')
    ? cancers
    : [{ label: 'All', value: '' }].concat(cancers);
  const selectedCancer =
    cancerTypes.find((c) => c.value === form.cancer) ?? cancerTypes[0];
  const leafColorOptions = isUser ? userColorOptions : colorOptions;

  function handleSearch(e) {
    onChange({ searchSamples: e });
  }

  function filterSampleOptions(inputValue = '', limit = 100) {
    return (attributes || [])
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
        {!isUser && (
          <Col md="auto">
            <Form.Group controlId="cancerType" className="mb-3">
              <Form.Label>Cancer Type</Form.Label>
              <Select
                name="cancerType"
                value={selectedCancer}
                options={cancerTypes}
                onChange={(e) => onChange({ cancer: e.value })}
                aria-label="Cancer Type Selector"
              />
            </Form.Group>
          </Col>
        )}
        <Col md="auto">
          <Form.Group controlId="color" className="mb-3">
            <Form.Label>Leaf Property</Form.Label>
            <Select
              name="color"
              value={form.color}
              options={leafColorOptions}
              onChange={(e) => onChange({ color: e })}
              aria-label="Leaf Property Selector"
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
              aria-label="Sample Search Selector"
            />
          </Form.Group>
        </Col>
      </Row>
    </Form>
  );
}
