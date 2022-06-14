import { useRecoilState, useRecoilValue } from 'recoil';
import { formState } from './mutProfiles.state';
import { Form, Row, Col } from 'react-bootstrap';
import CustomSelect from '../../../controls/select/select-old';
import Select from 'react-select';

export default function TreeLeafForm() {
  const [form, setForm] = useRecoilState(formState);
  const mergeForm = (state) => setForm({ ...form, ...state });

  function handleChange(event) {
    let { name, value, checked, type } = event.target;
    mergeForm({
      [name]: type === 'checkbox' ? checked : value,
    });
  }

  const options = [
    {
      label: 'SBS96',
      value: {
        profile: 'SBS96',
        signature_set: 'COSMIC_v3_Signatures_GRCh37_SBS96',
        signature: 'SBS84',
      },
    },
    {
      label: 'DBS78',
      value: {
        profile: 'DBS78',
        signature_set: 'COSMIC_v3_Signatures_GRCh37_DBS78',
        signature: 'DBS1',
      },
    },
    {
      label: 'ID1',
      value: {
        profile: 'ID83',
        signature_set: 'COSMIC_v3_Signatures_GRCh37_ID83',
        signature: 'ID1',
      },
    },
  ];

  return (
    <Form>
      <Row>
        <Col md="auto">
          <Form.Group controlId="profile" className="mb-3">
            <Form.Label>Profile</Form.Label>
            <Select
              name="option"
              // defaultValue={profileOptions[0]}
              options={options}
              onChange={(e) => mergeForm({ option: e })}
            />
          </Form.Group>
        </Col>
      </Row>
    </Form>
    // <Form className="p-3">
    //   <Row>
    //     <Col lg="auto">
    //       <CustomSelect
    //         id="mpSampleName"
    //         label="Sample Name"
    //         options={nameOptions}
    //         onChange={handleChange}
    //       />
    //     </Col>
    //     <Col lg="auto">
    //       <CustomSelect
    //         id="mpProfileType"
    //         label="Profile Type"
    //         value={selectProfile}
    //         options={profileOptions}
    //         onChange={handleProfile}
    //       />
    //     </Col>
    //     <Col lg="auto">
    //       <CustomSelect
    //         id="mpMatrixSize"
    //         label="Matrix Size"
    //         value={selectMatrix}
    //         options={matrixOptions}
    //         onChange={handleMatrix}
    //       />
    //     </Col>
    //     <Col lg="auto">
    //       <CustomSelect
    //         id="mpFilter"
    //         label="Filter"
    //         value={selectFilter}
    //         options={filterOptions}
    //         onChange={handleTag}
    //         disabled={source == 'public'}
    //       />
    //     </Col>
    //   </Row>
    // </Form>
  );
}
