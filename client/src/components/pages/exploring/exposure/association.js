import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchExpAssociation,
} from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';
import Select from '../../../controls/select/select';

const { Group, Label, Control, Check, Text } = Form;

export default function Association({ calculateAssociation }) {
  const { signatureNameOptions, userNameOptions, source } = useSelector(
    (state) => state.expExposure
  );
  const {
    toggleCancer,
    both,
    signatureName1,
    signatureName2,
    plotPath,
    plotURL,
    debugR,
    err,
    loading,
  } = useSelector((state) => state.expAssociation);
  const { projectID } = useSelector((state) => state.visualizeResults);

  useEffect(() => {
    if (plotPath) setRPlot(plotPath);
    else clearPlot();
  }, [plotPath]);

  useEffect(() => {
    if (source == 'public') {
      dispatchExpAssociation({
        signatureName1: signatureNameOptions[0],
        signatureName2: signatureNameOptions[1],
      });
    } else {
      dispatchExpAssociation({
        signatureName1: userNameOptions[0],
        signatureName2: userNameOptions[1],
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
          dispatchExpAssociation({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchExpAssociation({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) URL.revokeObjectURL(plotURL);
    dispatchExpAssociation({ plotPath: '', plotURL: '' });
  }

  return (
    <div>
      <Form>
        <LoadingOverlay active={loading} />
        <Row className="">
          <Col sm="6" md="3" className="my-auto">
            <Group controlId="toggleCancerType" className="d-flex">
              <Check
                type="checkbox"
                label="Selected Cancer Type Only"
                value={toggleCancer}
                checked={toggleCancer}
                onChange={() =>
                  dispatchExpAssociation({ toggleCancer: !toggleCancer })
                }
              />
            </Group>
          </Col>
          <Col sm="6" md="4" className="my-auto">
            <Group controlId="toggleBothSamples" className="d-flex">
              <Check
                type="checkbox"
                label="Number of Mutations assigned to both signature > 0"
                value={both}
                checked={both}
                onChange={(e) => dispatchExpAssociation({ both: !both })}
              />
            </Group>
          </Col>
          <Col sm="4" md="2">
            <Select
              className="mb-0"
              id="associationSignatureName1"
              label="Signature Name 1"
              value={signatureName1}
              options={
                source == 'public' ? signatureNameOptions : userNameOptions
              }
              onChange={(name) =>
                dispatchExpAssociation({ signatureName1: name })
              }
            />
          </Col>
          <Col sm="4" md="2">
            <Select
              className="mb-0"
              id="associationSignatureName2"
              label="Signature Name 2"
              value={signatureName2}
              options={
                source == 'public' ? signatureNameOptions : userNameOptions
              }
              onChange={(name) =>
                dispatchExpAssociation({ signatureName2: name })
              }
            />
          </Col>
          <Col sm="4" md="1" className="d-flex justify-content-end mt-auto">
            <Button variant="primary" onClick={calculateAssociation}>
              Calculate
            </Button>
          </Col>
        </Row>
        <div id="exposureAssociationPlot">
          {err && (
            <div>
              <p>
                An error has occured. Check the debug section for more info.
              </p>
              <p>Error: {err}</p>
            </div>
          )}
          {plotURL && (
            <Plot
              plotName={plotPath.split('/').slice(-1)[0]}
              plotURL={plotURL}
            />
          )}
          <Debug msg={debugR} />
        </div>
      </Form>
    </div>
  );
}
