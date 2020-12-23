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

const { Group, Check } = Form;

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
  const { projectID } = useSelector((state) => state.expExposure);

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
          <Col lg="3">
            <Group controlId="toggleCancerType">
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
          <Col lg="3">
            <Group controlId="toggleBothSamples">
              <Check
                type="checkbox"
                label="Number of Mutations assigned to both signature > 0"
                value={both}
                checked={both}
                onChange={(e) => dispatchExpAssociation({ both: !both })}
              />
            </Group>
          </Col>
          <Col lg="2">
            <Select
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
          <Col lg="2">
            <Select
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
          <Col lg="2" className="d-flex justify-content-end">
            <Button
              className="mt-auto mb-3"
              variant="primary"
              onClick={calculateAssociation}
            >
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
