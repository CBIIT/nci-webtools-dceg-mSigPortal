import { useRecoilState, useRecoilValue } from 'recoil';
import { useSelector } from 'react-redux';
import Select from 'react-select';
import { Form, Row, Col } from 'react-bootstrap';
import MultiSelect from '../../../controls/select/multiSelect';
import { formState, graphDataSelector, colorOptions, defaultFormState } from './treeLeaf.state';
import { useEffect } from 'react';

export default function TreeLeafForm() {
  const store = useSelector((state) => state.visualization);
  const study = store?.publicForm?.study?.value || 'PCAWG';
  const strategy = store?.publicForm?.strategy?.value || 'WGS';
  const cancers = store?.publicForm?.cancers?.filter(c => c.value !== '*ALL') || [];
  const cancerTypes = [{ label: 'All', value: '' }].concat(cancers);
  const [form, setForm] = useRecoilState(formState);
  const signatureSetName = form?.signatureSetName;
  const profile = form?.profile;
  const matrix = form?.matrix;
  const defaultCancer = cancerTypes.find(c => c.value === store?.publicForm?.cancer.value) || cancerTypes[0];
  const cancer = form?.cancerType?.value;
  const params = { study, strategy,  signatureSetName, profile, matrix, cancer };
  const { attributes } = useRecoilValue(graphDataSelector(params));
  const mergeForm = (state) => setForm({ ...form, ...state });
  useEffect(() => setForm({...defaultFormState, cancerType: defaultCancer}), [setForm, study, defaultCancer]); // reset form when study changes

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
          <Form.Group controlId="cancerType" className="mb-3">
            <Form.Label>Cancer Type</Form.Label>
            <Select
              name="cancerType"
              defaultValue={defaultCancer}
              options={cancerTypes}
              onChange={(e) => mergeForm({ cancerType: e })}
            />
          </Form.Group>
        </Col>        
        <Col md="auto">
          <Form.Group controlId="color" className="mb-3">
            <Form.Label>Leaf Property</Form.Label>
            <Select
              name="color"
              defaultValue={colorOptions[0]}
              options={colorOptions}
              onChange={(e) => mergeForm({ color: e })}
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
