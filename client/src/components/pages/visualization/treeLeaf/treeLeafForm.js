import { useRecoilState, useRecoilValue } from 'recoil';
import { useSelector } from 'react-redux';
import Select from 'react-select';
import { Form, Row, Col } from 'react-bootstrap';
import MultiSelect from '../../../controls/select/multiSelect';
import { formState, graphDataSelector, colorOptions } from './treeLeaf.state';

export default function TreeLeafForm() {
  const store = useSelector((state) => state.visualization);
  const study = store?.publicForm?.study?.value || 'PCAWG';
  const strategy = store?.publicForm?.strategy?.value || 'WGS';
  const signatureSetName = 'COSMIC_v3_Signatures_GRCh37_SBS96' 
  const profileMatrix = ["SBS96", "DBS78", "ID83"];
  const params = { study, strategy, signatureSetName, profileMatrix };
  const { attributes } = useRecoilValue(graphDataSelector(params));
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

  const customStyles = {
    control: (base, state) => ({
      ...base,
      background: '#f1e4ef',
      // match with the menu
      borderRadius: state.isFocused ? '3px 3px 0 0' : 3,
      // Overwrittes the different states of border
      borderColor: state.isFocused ? '#f1e4ef' : '#8e4b86',
      // Removes weird border around container
      boxShadow: state.isFocused ? null : null,
      '&:hover': {
        // Overwrittes the different states of border
        borderColor: state.isFocused ? '#8e4b86' : '#f1e4ef',
      },
    }),
    menu: (base) => ({
      ...base,
      // override border radius to match the box
      borderRadius: 0,
      // kill the gap
      marginTop: 0,
    }),
    menuList: (base) => ({
      ...base,
      // kill the white space on first and last option
      padding: 0,
    }),
  };

  return (
    <Form>
      <Row>
        <Col md="auto">
          <Form.Group controlId="color" className="mb-3">
            <Form.Label>Color By</Form.Label>
            <Select
              name="color"
              defaultValue={colorOptions[0]}
              options={colorOptions}
              onChange={(e) => mergeForm({ color: e })}
              styles={customStyles}
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
              styles={customStyles}
            />
          </Form.Group>
        </Col>
      </Row>
    </Form>
  );
}
