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

export function CosineSimilarity({ setPlot, submitR }) {
  const { pyTab, csPlotURL, cosineSimilarity } = useSelector(
    (state) => state.visualizeResults
  );

  const { nameOptions, profileOptions, matrixOptions } = pyTab;

  const {
    profileType1,
    matrixSize,
    profileType2,
    signatureSet,
    signatureSetOptions,
    selectName2,
    selectSigFormula,
    sigFormula,
    rPlots,
    rPlotIndex,
    submitOverlay,
    refSigOverlay,
    debugR,
  } = cosineSimilarity;

  // load first r plot after they are recieved
  useEffect(() => {
    if (rPlots.length && !rPlotIndex.length) {
      setPlot(0, 'cosineSimilarity');
    }
  }, [rPlots]);

  useEffect(() => {
    getSignatureSet(profileType2);
  }, [profileType2]);

  // get Signature Reference Sets
  async function getSignatureSet(profileType) {
    dispatchVisualizeResults({
      cosineSimilarity: { ...cosineSimilarity, refSigOverlay: true },
    });
    try {
      const response = await fetch(
        `${root}api/visualizeR/getSignatureReferenceSets`,
        {
          method: 'POST',
          headers: {
            Accept: 'application/json',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({ profileType: profileType }),
        }
      );
      const signatureSetOptions = await response.json();

      dispatchVisualizeResults({
        cosineSimilarity: {
          ...cosineSimilarity,
          signatureSetOptions: signatureSetOptions,
          signatureSet: signatureSetOptions[0],
          refSigOverlay: false,
        },
      });
    } catch (err) {
      dispatchError(err);
    } finally {
      //   dispatchVisualizeResults({ cosineSimilarity: { ...cosineSimilarity, refSigOverlay: false } });
    }
  }

  return (
    <div>
      <Form>
        <Label>Cosine Similarity Within Samples</Label>
        <div className="border rounded p-2">
          <Row className="justify-content-center">
            <Col sm="6">
              <Group controlId="profileType1">
                <Label>Profile Type</Label>
                <Control
                  as="select"
                  value={profileType1}
                  onChange={(e) =>
                    dispatchVisualizeResults({
                      cosineSimilarity: {
                        ...cosineSimilarity,
                        profileType1: e.target.value,
                      },
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
            <Col sm="6">
              <Label>Matrix Size</Label>
              <Control
                as="select"
                value={matrixSize}
                onChange={(e) =>
                  dispatchVisualizeResults({
                    cosineSimilarity: {
                      ...cosineSimilarity,
                      matrixSize: e.target.value,
                    },
                  })
                }
                custom
              >
                <option value="0" disabled>
                  Select
                </option>
                {matrixOptions.map((matrix, index) => {
                  return (
                    <option key={index} value={matrix}>
                      {matrix}
                    </option>
                  );
                })}
              </Control>
            </Col>
          </Row>
        </div>
      </Form>

      <Form>
        <Label>Cosine Similarity to Reference Signatures</Label>
        <LoadingOverlay active={refSigOverlay} />
        <div className="border rounded p-2">
          <Row className="justify-content-center">
            <Col sm="6">
              <Group controlId="profileType2">
                <Label>Profile Type</Label>
                <Control
                  as="select"
                  value={profileType2}
                  onChange={(e) => {
                    dispatchVisualizeResults({
                      cosineSimilarity: {
                        ...cosineSimilarity,
                        profileType2: e.target.value,
                      },
                    });
                  }}
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
            <Col sm="6">
              <Group controlId="signatureSet">
                <Label>Reference Signature Set</Label>
                <Control
                  disabled={!signatureSetOptions.length}
                  as="select"
                  value={signatureSet}
                  onChange={(e) =>
                    dispatchVisualizeResults({
                      cosineSimilarity: {
                        ...cosineSimilarity,
                        signatureSet: e.target.value,
                      },
                    })
                  }
                  custom
                >
                  <option value="0">Select</option>
                  {signatureSetOptions.map((signatureSet, index) => {
                    return (
                      <option key={index} value={signatureSet}>
                        {signatureSet}
                      </option>
                    );
                  })}
                </Control>
              </Group>
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

      {/* <Form>
        <Label>Additional Plots</Label>
        <div className="border rounded p-2">
          <LoadingOverlay active={submitOverlay} />
          <Row className="justify-content-center">
            <Col sm="4">
              <Group controlId="selectName2">
                <Label>Sample Name</Label>
                <Control
                  as="select"
                  value={selectName2}
                  onChange={(e) =>
                    dispatchVisualizeResults({
                      cosineSimilarity: {
                        ...cosineSimilarity,
                        selectName2: e.target.value,
                      },
                    })
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
                      cosineSimilarity: {
                        ...cosineSimilarity,
                        selectSigFormula: e.target.value,
                      },
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
                  dispatchVisualizeResults({
                    cosineSimilarity: {
                      ...cosineSimilarity,
                      sigFormula: e.target.value,
                    },
                  })
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
      </Form> */}

      <div
        className="mt-2 p-2 border rounded"
        style={{ display: rPlots.length ? 'block' : 'none' }}
      >
        <Col sm="auto" className="p-0">
          <Control
            as="select"
            value={rPlotIndex}
            onChange={(e) => setPlot(e.target.value, 'cosineSimilarity')}
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
            <img className="w-100 my-4" src={csPlotURL}></img>
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
