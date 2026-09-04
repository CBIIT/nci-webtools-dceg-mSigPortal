import { useRecoilState, useRecoilValue } from 'recoil';
import { useSelector } from 'react-redux';
import Select from 'react-select';
import { Form, Row, Col } from 'react-bootstrap';
import MultiSelect from '../../../controls/select/multiSelect';
import {
  formState,
  graphDataSelector,
  colorOptions,
  userColorOptions,
  defaultFormState,
  defaultUserFormState,
} from './treeLeaf.state';
import { useEffect } from 'react';

export default function TreeLeafForm({ state = {} }) {
  const { id, source } = state;
  const isUser = source === 'user';

  const publicForm = useSelector((store) => store.visualization.publicForm);
  const study = publicForm?.study?.value || 'PCAWG';
  const strategy = publicForm?.strategy?.value || 'WGS';
  const cancerValue = publicForm?.cancer?.value ?? '';
  const cancers = publicForm?.cancers?.filter((c) => c.value !== '*ALL') || [];
  const cancerTypes = [{ label: 'All', value: '' }].concat(cancers);
  const [form, setForm] = useRecoilState(formState);
  const signatureSetName = form?.signatureSetName;
  const profile = form?.profile || 'SBS';
  const matrix = form?.matrix || 96;
  const defaultCancer =
    cancerTypes.find((c) => c.value === cancerValue) || cancerTypes[0];
  const cancer = form?.cancerType?.value;

  const params = isUser
    ? { userId: id, source: 'user', profile: 'SBS', matrix: 96 }
    : { study, strategy, signatureSetName, profile, matrix, cancer };

  const graphData = useRecoilValue(graphDataSelector(params));
  if (graphData?.error) {
    throw new Error(graphData.error);
  }
  if (graphData === null && !(isUser && !id)) {
    throw new Error(
      isUser
        ? 'Failed to load Tree and Leaf data for this user session.'
        : 'Failed to load Tree and Leaf data for the selected study.'
    );
  }
  const { attributes } = graphData || {};
  const mergeForm = (next) => setForm({ ...form, ...next });
  const leafColorOptions = isUser ? userColorOptions : colorOptions;

  // Reset form when source/study/session changes. Depend on primitive values only —
  // object identities like defaultCancer change every render and would loop setForm.
  useEffect(() => {
    if (isUser) {
      setForm({ ...defaultUserFormState });
      return;
    }
    const cancerOption = [{ label: 'All', value: '' }]
      .concat((publicForm?.cancers || []).filter((c) => c.value !== '*ALL'))
      .find((c) => c.value === cancerValue) || { label: 'All', value: '' };
    setForm({ ...defaultFormState, cancerType: cancerOption });
  }, [setForm, study, cancerValue, isUser, id]);

  function handleSearch(e) {
    mergeForm({ searchSamples: e });
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
                value={form.cancerType || defaultCancer}
                options={cancerTypes}
                onChange={(e) => mergeForm({ cancerType: e })}
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
              onChange={(e) => mergeForm({ color: e })}
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
