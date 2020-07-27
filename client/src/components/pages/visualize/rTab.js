import React, { useEffect, useState } from 'react';
import { Form, Row, Col, Button, Card, Nav } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchVisualizeResults,
} from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

const { Group, Label, Control } = Form;

const root =
  process.env.NODE_ENV === 'development'
    ? 'http://localhost:8330/'
    : window.location.pathname;

export default function RTab({ setPlot, submitR }) {
  const { pyTab, rPlotURL, rTab } = useSelector(
    (state) => state.visualizeResults
  );

  const { nameOptions, profileOptions } = pyTab;

  const {
    sigProfileType,
    matrixSize,
    signatureSet,
    selectName2,
    selectSigFormula,
    sigFormula,
    rPlots,
    rPlotIndex,
    submitOverlay,
    debugR,
  } = rTab;

  // load first r plot after they are recieved
  useEffect(() => {
    if (rPlots.length && !rPlotIndex.length) {
      setPlot(0, 'r');
    }
  }, [rPlots]);

  return (
    <div>
      <Form>
        <Label>Additional Plots</Label>
        <div className="border rounded p-2">
          <LoadingOverlay active={submitOverlay} />
          <Row className="justify-content-center">
            <Col sm="4">
              <Group controlId="sigProfileType">
                <Label>Signature Profile Type</Label>
                <Control
                  as="select"
                  value={sigProfileType}
                  onChange={(e) =>
                    dispatchVisualizeResults({
                      sigProfileType: e.target.value,
                    })
                  }
                  custom
                >
                  {profileOptions.map((profile, index) => {
                    return (
                      <option key={index} value={profile}>
                        {profile}
                      </option>
                    );
                  })}
                </Control>
              </Group>
            </Col>
            <Col sm="4">
              <Group controlId="signatureSet">
                <Label>Reference Set</Label>
                <Control
                  as="select"
                  value={signatureSet}
                  onChange={(e) =>
                    dispatchVisualizeResults({ signatureSet: e.target.value })
                  }
                  custom
                >
                  <option value="COSMIC v3 Signatures (SBS)">
                    COSMIC v3 Signatures (SBS)
                  </option>
                </Control>
              </Group>
            </Col>
            <Col sm="4">
              <Group controlId="selectName2">
                <Label>Sample Name</Label>
                <Control
                  as="select"
                  value={selectName2}
                  onChange={(e) =>
                    dispatchVisualizeResults({ selectName2: e.target.value })
                  }
                  custom
                >
                  <option value="0" disabled>
                    Select
                  </option>
                  {nameOptions.map((sampleName, index) => {
                    return (
                      <option key={index} value={sampleName}>
                        {sampleName}
                      </option>
                    );
                  })}
                </Control>
              </Group>
            </Col>
          </Row>
          <Row className="justify-content-center">
            <Col sm="6">
              <Group
                controlId="selectSigFormula"
                className="d-flex align-items-center"
              >
                <Label className="mr-auto">Signature/Formula</Label>
                <Control
                  as="select"
                  value={selectSigFormula}
                  onChange={(e) =>
                    dispatchVisualizeResults({
                      selectSigFormula: e.target.value,
                    })
                  }
                  custom
                  style={{ width: '150px' }}
                >
                  <option value="signature">Signature</option>
                  <option value="formula">Formula</option>
                </Control>
              </Group>
            </Col>
            <Col sm="6">
              <Control
                type="text"
                placeholder={selectSigFormula}
                value={sigFormula}
                onChange={(e) =>
                  dispatchVisualizeResults({ sigFormula: e.target.value })
                }
              ></Control>
            </Col>
          </Row>
          <Row>
            <Col>
              <Button variant="primary" onClick={() => submitR()}>
                Calculate
              </Button>
            </Col>
          </Row>
        </div>
      </Form>
      <div
        className="mt-2 p-2 border rounded"
        style={{ display: rPlots.length ? 'block' : 'none' }}
      >
        <Col sm="auto" className="p-0">
          <Control
            as="select"
            value={rPlotIndex}
            onChange={(e) => setPlot(e.target.value, 'r')}
            custom
          >
            {rPlots.map((plot, index) => {
              return (
                <option key={index} value={index}>
                  {plot}
                </option>
              );
            })}
          </Control>
        </Col>
        <Row>
          <Col>
            <img className="w-100 my-4" src={rPlotURL}></img>
          </Col>
        </Row>
      </div>
      <div className="border rounded p-1 my-2">
        <div>R output</div>
        <div className="border">
          {Array.isArray(debugR) ? (
            debugR.map((line, index) => {
              return (
                <p key={index} className="m-0">
                  [{index}] {line}
                </p>
              );
            })
          ) : (
            <p>{debugR}</p>
          )}
        </div>
      </div>
    </div>
  );
}
