import { useRecoilState, useRecoilValue } from "recoil";
import { formState } from "./mutProfiles.state";
import { Form, Row, Col } from "react-bootstrap";
import CustomSelect from "../../../controls/select/select";
import Select from "react-select";

export default function TreeLeafForm() {
  const [form, setForm] = useRecoilState(formState);
  const mergeForm = (state) => setForm({ ...form, ...state });

  function handleChange(event) {
    let { name, value, checked, type } = event.target;
    mergeForm({
      [name]: type === "checkbox" ? checked : value,
    });
  }

  const profileOptions = [
    {
      label: "SBS96",
      value: {
        Profile: "SBS96",
        Signature_set_name: "COSMIC_v3_Signatures_GRCh37_SBS96",
        Signature_name: "SBS84",
      },
    },
    {
      label: "DBS78",
      value: {
        Profile: "DBS78",
        Signature_set_name: "COSMIC_v3_Signatures_GRCh37_DBS78",
        Signature_name: "DBS1",
      },
    },
    {
      label: "ID1",
      value: {
        Profile: "ID83",
        Signature_set_name: "COSMIC_v3_Signatures_GRCh37_ID83",
        Signature_name: "ID1",
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
              name="profile"
              // defaultValue={profileOptions[0]}
              options={profileOptions}
              onChange={(e) => mergeForm({ profile: e })}
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
