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

  useEffect(() => {
    if (plotPath) setRPlot(plotPath);
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
        const response = await fetch(`api/results/${plotPath}`);
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

  return (
    <div>
      <Form>
        <LoadingOverlay active={loading} />
        <div className="px-4">
          <Row className="justify-content-center">
            <Col sm="2" className="my-auto">
              <Group controlId="toggleCancerType" className="d-flex">
                <Label className="mr-4">Use Cancer Type</Label>
                <Check inline id="toggleCancerType">
                  <Check.Input
                    type="checkbox"
                    value={toggleCancer}
                    checked={toggleCancer}
                    onChange={() =>
                      dispatchExpAssociation({ toggleCancer: !toggleCancer })
                    }
                  />
                </Check>
              </Group>
            </Col>
            <Col sm="2" className="my-auto">
              <Group controlId="toggleBothSamples" className="d-flex">
                <Label className="mr-4">Both Signatures</Label>
                <Check inline id="toggleBothSamples">
                  <Check.Input
                    type="checkbox"
                    value={both}
                    checked={both}
                    onChange={(e) => dispatchExpAssociation({ both: !both })}
                  />
                </Check>
              </Group>
            </Col>
            <Col sm="2">
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
            <Col sm="2">
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
            <Col sm="3" />
            <Col sm="1" className="my-auto">
              <Button variant="primary" onClick={calculateAssociation}>
                Calculate
              </Button>
            </Col>
          </Row>
          <div id="associationPlot">
            {err && (
              <p>
                An error has occured. Check the debug section for more info.
              </p>
            )}
            {plotURL && (
              <Plot
                plotName={plotPath.split('/').slice(-1)[0]}
                plotURL={plotURL}
              />
            )}
            <Debug msg={debugR} />
          </div>
        </div>
      </Form>
    </div>
  );
}
