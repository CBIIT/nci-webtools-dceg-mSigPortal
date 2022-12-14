import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import Select from '../../../controls/select/selectForm';
import Description from '../../../controls/description/description';
import { useSelector, useDispatch } from 'react-redux';
import { actions as exposureActions } from '../../../../services/store/exposure';

import { NavHashLink } from 'react-router-hash-link';

const actions = { ...exposureActions };

export default function MsBurdenForm() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.exposure);

  const mergeMsBurden = (state) =>
    dispatch(actions.mergeExposure({ msBurden: state }));

  const { signatureNames } = store.main;
  const { signatureName } = store.msBurden;

  const signatureNameOptions = signatureNames.length
    ? [...new Set(signatureNames.map((d) => d))]
        .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
        .map((e) => ({
          label: e,
          value: e,
        }))
    : [];

  const { control, setValue, watch } = useForm({
    defaultValues: signatureName,
  });

  // set inital
  useEffect(() => {
    if (!signatureName && signatureNameOptions.length) {
      mergeMsBurden({ signatureName: signatureNameOptions[0] });
    }
  }, [signatureNameOptions]);

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
    <div>
      <hr />
      <Form className="p-3">
        <Row>
          <Col lg="auto">
            <Select
              name="signatureName"
              label="Signature Name"
              value={signatureName}
              control={control}
              options={signatureNameOptions}
              onChange={(name) => mergeMsBurden({ signatureName: name })}
              styles={customStyles}
              //onChange={handleSignatureName}
            />
          </Col>
          {/* <Col lg="auto" className="d-flex">
            <Button
              className="mt-auto mb-2"
              variant="primary"
              //onClick={}
            >
              Recalculate
            </Button>
          </Col> */}
        </Row>
      </Form>
    </div>
  );
}
