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

export default function MsAssociationForm() {
  const dispatch = useDispatch();
  const mergeMsAssociation = (state) =>
    dispatch(actions.mergeExposure({ msAssociation: state }));
  const exposure = useSelector((state) => state.exposure);
  const { publicForm, main, msAssociation } = useSelector(
    (state) => state.exposure
  );
  const [signatureOptionParams, setSignatureOptionParams] = useState('');
  const [params, setParams] = useState(null);
  const {
    data: signatureNameOptions,
    error: signatureNameError,
    isFetching: signatureNameIsFetching,
  } = useMsAssociationOptionsQuery(signatureOptionParams, {
    skip: !signatureOptionParams,
  });

  const { data, error, isFetching } = useMsAssociationQuery(params, {
    skip: !params,
  });

  const { control, handleSubmit, watch } = useForm({
    defaultValues: msAssociation,
  });

  const { signatureName1, signatureName2 } = watch();

  //   console.log(signatureNameOptions);
  //   console.log(msAssociation);
  //   console.log(publicForm);
  // query signature name options
  useEffect(() => {
    if (publicForm.study) {
      setSignatureOptionParams({
        study: publicForm.study.value,
        strategy: publicForm.strategy.value,
        signatureSetName: publicForm.signatureSetName.value,
        ...(!publicForm.useAllCancer && { cancer: publicForm.cancer.value }),
      });
    } else if (main.id) {
      setSignatureOptionParams({ userId: main.id });
    }
  }, [publicForm, main]);

  function getSignatureNameOptions(name) {
    return name ? name.map((e) => ({ label: e, value: e })) : [];
  }

  function onSubmit(data) {
    console.log(data);
    mergeMsAssociation(data);
    const params = {
      study: publicForm.study.value,
      strategy: publicForm.strategy.value,
      signatureSetName: publicForm.signatureSetName.value,
      cancer: publicForm.cancer.value,

      signatureName:
        data.signatureName1.value + ';' + data.signatureName2.value,
    };
    setParams(params);
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
              value={msAssociation.signatureName1}
              disabled={!signatureNameOptions}
              control={control}
              options={signatureNameOptions}
              //onChange={(name) => mergeMsAssociation({ signatureName1: name })}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="signatureName2"
              label="Signature Name 2"
              value={msAssociation.signatureName2}
              disabled={!signatureNameOptions}
              control={control}
              options={signatureNameOptions}
              //onChange={(name) => mergeMsAssociation({ signatureName2: name })}
            />
          </Col>
          <Col lg="auto">
            <Group controlId="toggleBothSamples">
              <Check
                type="checkbox"
                label="Samples Detected Both Signatures"
                value={msAssociation.both}
                checked={msAssociation.both}
                // onChange={(e) =>
                //   mergeMsAssociation({ both: !msAssociation.both })
                // }
              />
            </Group>
          </Col>
          <Col lg="auto" className="d-flex">
            <Button
              className="mt-auto mb-3"
              //disabled={!signatureName1 || !signatureName2}
              variant="primary"
              type="submit"
            >
              Calculate
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
