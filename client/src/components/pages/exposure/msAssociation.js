import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { actions as exposureActions } from '../../../services/store/exposure';
import { actions as modalActions } from '../../../services/store/modal';
import Plot from '../../controls/plot/plot';
import Select from '../../controls/select/select';
import Debug from '../../controls/debug/debug';

const actions = { ...exposureActions, ...modalActions };
const { Group, Check } = Form;

export default function MsAssociation({ calculateAssociation }) {
  const dispatch = useDispatch();
  const exposure = useSelector((state) => state.exposure);
  const {
    both,
    signatureName1,
    signatureName2,
    plotPath,
    debugR,
    err,
    loading,
  } = exposure.msAssociation;
  const {
    projectID,
    signatureNameOptions,
    userNameOptions,
    source,
  } = exposure.exposureState;
  const mergeExposure = (state) =>
    dispatch(actions.mergeExposure({ exposure: state }));
  const mergeMsAssociation = (state) =>
    dispatch(actions.mergeExposure({ msAssociation: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  // apply default signature names
  useEffect(() => {
    if (source == 'public') {
      mergeMsAssociation({
        signatureName1: signatureNameOptions[0],
        signatureName2: signatureNameOptions[1],
      });
    } else {
      mergeMsAssociation({
        signatureName1: userNameOptions[0] || '',
        signatureName2: userNameOptions[1] || '',
      });
    }
  }, [signatureNameOptions, userNameOptions, source]);

  return (
    <div>
      <div className="p-3">
        <b>Mutational Signature Association</b>
        <p className="m-0">
          The scatter plot below illustrates the associations between two
          selected. On the x-axis is the number of mutations (log 10) assigned
          to Signature Name 1, and on the y-axis is the number of mutations
          (log10) assigned to Signature Name 2. Use the parameter on top panel
          “Number of Mutations assigned to both signature >0” to remove the
          samples without assigning any mutations to these selected two
          signatures. Use the parameter on the left panel “Cancer Type Only” to
          perform association analysis on selected the cancer type only.
        </p>
      </div>
      <hr />
      <Form className="p-3">
        <LoadingOverlay active={loading} />
        <Row>
          <Col lg="3">
            <Group controlId="toggleBothSamples">
              <Check
                type="checkbox"
                label="Number of Mutations assigned to both signature > 0"
                value={both}
                checked={both}
                onChange={(e) => mergeMsAssociation({ both: !both })}
              />
            </Group>
          </Col>
          <Col lg="2">
            <Select
              disabled={source == 'user' && !userNameOptions.length}
              id="associationSignatureName1"
              label="Signature Name 1"
              value={signatureName1}
              options={
                source == 'public' ? signatureNameOptions : userNameOptions
              }
              onChange={(name) => mergeMsAssociation({ signatureName1: name })}
            />
          </Col>
          <Col lg="2">
            <Select
              disabled={source == 'user' && !userNameOptions.length}
              id="associationSignatureName2"
              label="Signature Name 2"
              value={signatureName2}
              options={
                source == 'public' ? signatureNameOptions : userNameOptions
              }
              onChange={(name) => mergeMsAssociation({ signatureName2: name })}
            />
          </Col>
          <Col />
          <Col lg="2" className="d-flex">
            <Button
              className="ml-auto mb-auto"
              disabled={
                !signatureName1 ||
                !signatureName2 ||
                (source == 'user' && !projectID)
              }
              variant="primary"
              onClick={calculateAssociation}
            >
              Calculate
            </Button>
          </Col>
        </Row>
      </Form>
      <div id="exposureAssociationPlot">
        {err && (
          <div>
            <hr />
            <p className="p-3 text-danger">{err}</p>
          </div>
        )}
        {plotPath && (
          <>
            <hr />
            <Plot
              className="p-3"
              title="Mutational Signature Association"
              downloadName={plotPath.split('/').slice(-1)[0]}
              plotPath={`api/results/${projectID}${plotPath}`}
              maxHeight="900px"
            />
          </>
        )}
        {/* <Debug msg={debugR} /> */}
      </div>
    </div>
  );
}
