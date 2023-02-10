import React, { useEffect, useState } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import Select from '../../../controls/select/selectForm';
import { useSelector, useDispatch } from 'react-redux';
import { actions as exposureActions } from '../../../../services/store/exposure';
import {
  useMsAssociationOptionsQuery,
  useMsAssociationQuery,
} from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plotly from '../../../controls/plotly/plot/plot';

const actions = { ...exposureActions };
const { Group, Check } = Form;

export default function MsAssociationForm({ state, form, setForm }) {
  const [signatureOptionParams, setSignatureOptionParams] = useState('');
  const [params, setParams] = useState(null);

  const [bothCheck, setBothCheck] = useState(false);
  const {
    data: signatureNameOptions,
    error: signatureNameError,
    isFetching: signatureNameIsFetching,
  } = useMsAssociationOptionsQuery(signatureOptionParams, {
    skip: !signatureOptionParams,
  });

  const { signatureName1, signatureName2 } = form;

  const { study, strategy, signatureSetName, cancer, useAllCancer, id } = state;

  const { data, error, isFetching } = useMsAssociationQuery(params, {
    skip: !params,
  });

  const { control, handleSubmit, watch, setValue } = useForm({
    defaultValues: form,
  });

  // set inital
  useEffect(() => {
    if (!form.signatureName1 && signatureNameOptions) {
      setForm({ ...form, signatureName1: signatureNameOptions[0] });
    }
    if (!form.signatureName2 && signatureNameOptions) {
      setForm({ ...form, signatureName2: signatureNameOptions[1] });
    }
  }, [signatureNameOptions]);

  useEffect(() => {
    if (study) {
      setSignatureOptionParams({
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        ...(!useAllCancer && { cancer: cancer.value }),
      });
    } else if (id) {
      console.log(id);
      setSignatureOptionParams({ userId: id });
    }
  }, [state]);

  // useEffect(() => {
  //   if (signatureName1 && signatureName2 && id && signatureNameOptions.length) {
  //     setParams({
  //       signatureName: signatureName1.value + ';' + signatureName2.value,
  //       both: bothCheck,
  //       userId: id,
  //     });
  //   } else if (
  //     signatureName1 &&
  //     signatureName2 &&
  //     study &&
  //     signatureNameOptions.length
  //   ) {
  //     setParams({
  //       signatureName: signatureName1.value + ';' + signatureName2.value,
  //       both: bothCheck,
  //       study: study.value,
  //       strategy: strategy.value,
  //       signatureSetName: signatureSetName.value,
  //       ...(!useAllCancer && { cancer: cancer.value }),
  //     });
  //   }
  // }, [signatureName1, signatureName2, id]);

  function handleSignatureName1(e) {
    setForm({ ...form, signatureName1: e });
  }

  function handleSignatureName2(e) {
    setForm({ ...form, signatureName2: e });
  }
  function onSubmit() {
    if (signatureName1 && signatureName2 && id) {
      setParams({
        signatureName: signatureName1.value + ';' + signatureName2.value,
        both: bothCheck,
        userId: id,
      });
    } else if (signatureName1 && signatureName2 && study) {
      setParams({
        signatureName: signatureName1.value + ';' + signatureName2.value,
        both: bothCheck,
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        ...(!useAllCancer && { cancer: cancer.value }),
      });
    }
  }

  return (
    <div>
      <Form className="p-3" onSubmit={handleSubmit(onSubmit)}>
        <LoadingOverlay active={signatureNameIsFetching} />
        <Row>
          <Col lg="auto">
            <Select
              name="signatureName1"
              label="Signature Name 1"
              value={form.signatureName1}
              disabled={!signatureNameOptions}
              control={control}
              options={signatureNameOptions}
              onChange={handleSignatureName1}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="signatureName2"
              label="Signature Name 2"
              value={form.signatureName2}
              disabled={!signatureNameOptions}
              control={control}
              options={signatureNameOptions}
              onChange={handleSignatureName2}
            />
          </Col>
          <Col lg="auto">
            <Group controlId="toggleBothSamples">
              <Check
                name="signatureCheck"
                type="checkbox"
                label="Samples Detected Both Signatures"
                defaultChecked={bothCheck}
                onChange={(e) => {
                  setBothCheck(e.target.checked);
                }}
              />
            </Group>
          </Col>
          <Col lg="auto" className="d-flex">
            <Button
              className="mt-auto mb-3"
              disabled={!signatureName1 || !signatureName2}
              variant="primary"
              type="submit"
            >
              Recalculate
            </Button>
          </Col>
        </Row>
      </Form>
      <div id="exposureAssociationPlot">
        {signatureNameError && (
          <div>
            <hr />
            <p className="p-3 text-danger">{signatureNameError}</p>
          </div>
        )}
        {data && (
          <Plotly
            className="w-100"
            data={data.traces}
            layout={data.layout}
            config={data.config}
          />
        )}
        {/* <Debug msg={debugR} /> */}
      </div>
    </div>
  );
}
