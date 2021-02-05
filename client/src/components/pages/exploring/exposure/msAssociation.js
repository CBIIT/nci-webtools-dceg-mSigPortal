import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { actions as exploringActions } from '../../../../services/store/exploring';
import { actions as modalActions } from '../../../../services/store/modal';
import Plot from '../../../controls/plot/plot';
import Select from '../../../controls/select/select';
import Debug from '../../../controls/debug/debug';

const actions = { ...exploringActions, ...modalActions };
const { Group, Check } = Form;

export default function MsAssociation({ calculateAssociation, handleSet }) {
  const dispatch = useDispatch();
  const exploring = useSelector((state) => state.exploring);
  const {
    both,
    signatureName1,
    signatureName2,
    plotPath,
    plotURL,
    debugR,
    err,
    loading,
  } = exploring.msAssociation;
  const {
    projectID,
    signatureNameOptions,
    userNameOptions,
    source,
    study,
    strategy,
    refSignatureSet,
  } = exploring.exposure;
  const mergeExploring = (state) =>
    dispatch(actions.mergeExploring({ exploring: state }));
  const mergeMsAssociation = (state) =>
    dispatch(actions.mergeExploring({ msAssociation: state }));
  const mergeError = (state) => dispatch(actions.mergeModal({ error: state }));

  useEffect(() => {
    plotPath ? setRPlot(plotPath) : clearPlot();
  }, [plotPath, err, debugR, projectID]);

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

  async function setRPlot(plotPath) {
    if (plotPath) {
      try {
        const response = await fetch(`api/results/${projectID}${plotPath}`);
        if (!response.ok) {
          // console.log(await response.json());
        } else {
          const pic = await response.blob();
          const objectURL = URL.createObjectURL(pic);

          if (plotURL) URL.revokeObjectURL(plotURL);
          mergeMsAssociation({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        mergeError({ visible: true, message: err.message });
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      mergeMsAssociation({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      mergeMsAssociation({ plotURL: '' });
    }
  }

  return (
    <div>
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
          <Col lg="2" className="d-flex justify-content-end">
            <Button
              disabled={!signatureName1 || !signatureName2}
              className="mt-auto mb-3"
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
              plotName="Mutational Signature Association"
              plotURL={plotURL}
            />
          </>
        )}
        {/* <Debug msg={debugR} /> */}
      </div>
    </div>
  );
}
