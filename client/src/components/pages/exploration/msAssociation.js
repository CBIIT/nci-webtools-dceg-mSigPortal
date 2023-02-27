import React from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import Description from '../../controls/description/description';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { actions as explorationActions } from '../../../services/store/exploration';
import { actions as modalActions } from '../../../services/store/modal';
import SvgContainer from '../../controls/svgContainer/svgContainer';
import CustomSelect from '../../controls/select/select-old';

const actions = { ...explorationActions, ...modalActions };
const { Group, Check } = Form;

export default function MsAssociation({ calculateAssociation }) {
  const dispatch = useDispatch();
  const exploration = useSelector((state) => state.exploration);
  const {
    both,
    signatureName1,
    signatureName2,
    plotPath,
    debugR,
    err,
    loading,
  } = exploration.msAssociation;
  const { id, signatureNameOptions, userNameOptions, source } =
    exploration.main;

  const mergeMsAssociation = (state) =>
    dispatch(actions.mergeExploration({ msAssociation: state }));

  return (
    <div>
      <div className="p-3">
        <b>Mutational Signature Association</b>
        <Description
          less="The scatter plot below illustrates the associations between two selected mutational signatures."
          more="On the x-axis is the number of mutations (log 10) assigned to Signature Name 1, and on the y-axis is the number of mutations (log10) assigned to Signature Name 2. Use the parameter from the top panel, [Samples Detected Both Signatures] to remove the samples that do not carry both signatures before running the association analysis. Use the parameter in the left panel [Cancer Type Only] to perform association analyses on the selected cancer type only."
        />
      </div>
      <hr />
      <Form className="p-3">
        <LoadingOverlay active={loading} />
        <Row>
          <Col lg="auto">
            <CustomSelect
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
          <Col lg="auto">
            <CustomSelect
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
          <Col lg="auto">
            <Group controlId="toggleBothSamples">
              <Check
                type="checkbox"
                label="Samples Detected Both Signatures"
                value={both}
                checked={both}
                onChange={(e) => mergeMsAssociation({ both: !both })}
              />
            </Group>
          </Col>
          <Col lg="auto" className="d-flex">
            <Button
              className="mt-auto mb-3"
              disabled={
                !signatureName1 ||
                !signatureName2 ||
                (source == 'user' && !id)
              }
              variant="primary"
              onClick={calculateAssociation}
            >
              Recalculate
            </Button>
          </Col>
        </Row>
      </Form>
      <div id="explorationAssociationPlot">
        {err && (
          <div>
            <hr />
            <p className="p-3 text-danger">{err}</p>
          </div>
        )}
        {plotPath && (
          <>
            <hr />
            <SvgContainer
              className="p-3"
              title="Mutational Signature Association"
              downloadName={plotPath.split('/').slice(-1)[0]}
              plotPath={`web/data/${plotPath}`}
              height="1100px"
            />
          </>
        )}
        {/* <Debug msg={debugR} /> */}
      </div>
    </div>
  );
}
