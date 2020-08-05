import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchVisualizeResults,
  dispatchProfileComparison,
} from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

const { Group, Label, Control } = Form;

export default function ProfileComparison({ submitR }) {
  const { nameOptions, profileOptions } = useSelector((state) => state.pyTab);
  const rootURL = window.location.pathname;
  const {
    withinProfileType,
    withinSampleName1,
    withinSampleName2,
    refProfileType,
    refSampleName,
    refSignatureSet,
    refSignatureSetOptions,
    refCompare,
    withinPlotPath,
    refPlotPath,
    withinPlotURL,
    refPlotURL,
    displayWithin,
    displayRefSig,
    debugR,
    submitOverlay,
  } = useSelector((state) => state.profileComparison);

  // load r plots after they are recieved
  useEffect(() => {
    if (withinPlotPath && withinPlotPath.length) {
      setRPlot(withinPlotPath, 'within');
    }
  }, [withinPlotPath]);

  useEffect(() => {
    if (refPlotPath && refPlotPath.length) {
      setRPlot(refPlotPath, 'refsig');
    }
  }, [refPlotPath]);

  // call r wrapper on load
  // useEffect(() => {
  //   if (!withinPlotURL.length && !refPlotURL.length) {
  //     calculateR('cosineSimilarityWithin', {
  //       profileType: profileType1,
  //       matrixSize: matrixSize.replace('-', ''),
  //     });
  //     calculateR('cosineSimilarityRefSig', {
  //       profileType: profileType2,
  //       signatureSet: signatureSet,
  //     });
  //   }
  // }, [withinPlotPath, refPlotPath]);

  // retrieve signature set options on change
  useEffect(() => {
    getSignatureSet(refProfileType);
  }, [refProfileType]);

  async function setRPlot(plotPath, type) {
    try {
      const response = await fetch(`${rootURL}visualize/svg`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ path: plotPath }),
      });
      if (!response.ok) {
        const { msg } = await response.json();
        dispatchError(msg);
      } else {
        const pic = await response.blob();
        const objectURL = URL.createObjectURL(pic);

        if (type == 'within') {
          if (withinPlotURL.length) URL.revokeObjectURL(withinPlotURL);
          dispatchProfileComparison({
            withinPlotURL: objectURL,
          });
        } else {
          if (refPlotURL.length) URL.revokeObjectURL(refPlotURL);
          dispatchProfileComparison({
            refPlotURL: objectURL,
          });
        }
      }
    } catch (err) {
      dispatchError(err);
    }
  }

  // get Signature Reference Sets for dropdown options
  async function getSignatureSet(profileType) {
    if (profileType.length) {
      dispatchProfileComparison({ submitOverlay: true });
      try {
        const response = await fetch(
          `${rootURL}api/visualizeR/getSignatureReferenceSets`,
          {
            method: 'POST',
            headers: {
              Accept: 'application/json',
              'Content-Type': 'application/json',
            },
            body: JSON.stringify({ profileType: profileType }),
          }
        );
        if (response.ok) {
          const signatureSetOptions = await response.json();

          dispatchProfileComparison({
            refSignatureSetOptions: signatureSetOptions,
            refSignatureSet: signatureSetOptions[0],
            submitOverlay: false,
          });
        } else {
          dispatchError(await response.json());
          dispatchProfileComparison({ submitOverlay: false });
        }
      } catch (err) {
        dispatchError(err);
        dispatchProfileComparison({ submitOverlay: false });
      }
    }
  }

  async function calculateR(fn, args) {
    dispatchProfileComparison({ submitOverlay: true, debugR: '' });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();
        dispatchProfileComparison({
          debugR: err,
          submitOverlay: false,
        });
      } else {
        const { debugR, output } = await response.json();

        if (fn == 'profileComparisonWithin') {
          dispatchProfileComparison({ withinPlotPath: '' });
          dispatchProfileComparison({
            debugR: debugR,
            submitOverlay: false,
            withinPlotPath: output.plotPath,
          });
        } else {
          {
            dispatchProfileComparison({ refPlotPath: '' });
            dispatchProfileComparison({
              debugR: debugR,
              submitOverlay: false,
              refPlotPath: output.plotPath,
            });
          }
        }
      }
    } catch (err) {
      dispatchError(err);
      dispatchProfileComparison({ submitOverlay: false });
    }
  }

  return (
    <div>
      <LoadingOverlay active={submitOverlay} />
      <Form>
        <Label>
          <Button
            variant="link"
            className="p-0 font-weight-bold"
            onClick={() =>
              dispatchProfileComparison({
                displayWithin: !displayWithin,
              })
            }
          >
            Cosine Similarity Within Samples
          </Button>
        </Label>
        <div
          className="border rounded p-2"
          style={{ display: displayWithin ? 'block' : 'none' }}
        >
          <Row className="justify-content-center">
            <Col sm="2">
              <Group controlId="profileTypeWithin">
                <Label>Profile Type</Label>
                <Control
                  as="select"
                  value={withinProfileType}
                  onChange={(e) =>
                    dispatchProfileComparison({
                      withinProfileType: e.target.value,
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
              <Label>Sample Name 1</Label>
              <Control
                as="select"
                value={withinSampleName1}
                onChange={(e) =>
                  dispatchProfileComparison({
                    withinSampleName1: e.target.value,
                  })
                }
                custom
              >
                <option value="0" disabled>
                  Select
                </option>
                {nameOptions.map((name, index) => {
                  return (
                    <option key={index} value={name}>
                      {name}
                    </option>
                  );
                })}
              </Control>
            </Col>
            <Col sm="4">
              <Label>Sample Name 2</Label>
              <Control
                as="select"
                value={withinSampleName2}
                onChange={(e) =>
                  dispatchProfileComparison({
                    withinSampleName2: e.target.value,
                  })
                }
                custom
              >
                <option value="0" disabled>
                  Select
                </option>
                {nameOptions.map((name, index) => {
                  return (
                    <option key={index} value={name}>
                      {name}
                    </option>
                  );
                })}
              </Control>
            </Col>
            <Col sm="2" className="m-auto">
              <Button
                variant="primary"
                onClick={() =>
                  calculateR('profileComparisonWithin', {
                    profileType: withinProfileType,
                    sampleName1: withinSampleName1,
                    sampleName2: withinSampleName2,
                  })
                }
              >
                Calculate
              </Button>
            </Col>
          </Row>

          <div
            id="pcWithinPlot"
            style={{ display: withinPlotURL.length ? 'block' : 'none' }}
          >
            <div className="d-flex">
              <a
                className="px-2 py-1"
                href={withinPlotURL}
                download={withinPlotURL.split('/').slice(-1)[0]}
              >
                Download Plot
              </a>
            </div>
            <div className="p-2 border rounded">
              <Row>
                <Col>
                  <img className="w-100 my-4" src={withinPlotURL}></img>
                </Col>
              </Row>
            </div>
          </div>
        </div>
      </Form>

      <Form className="my-2">
        <Label>
          <Button
            variant="link"
            className="p-0 font-weight-bold"
            onClick={() =>
              dispatchProfileComparison({
                displayRefSig: !displayRefSig,
              })
            }
          >
            Comparison to Reference Signatures
          </Button>
        </Label>
        <div
          className="border rounded p-2"
          style={{ display: displayRefSig ? 'block' : 'none' }}
        >
          <Row className="justify-content-center">
            <Col sm="2">
              <Group controlId="profileTypeRefSig">
                <Label>Profile Type</Label>
                <Control
                  as="select"
                  value={refProfileType}
                  onChange={(e) => {
                    dispatchProfileComparison({
                      refProfileType: e.target.value,
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
            <Col sm="2">
              <Group controlId="sampleNameRefSig">
                <Label>Sample Name</Label>
                <Control
                  as="select"
                  value={refSampleName}
                  onChange={(e) => {
                    dispatchProfileComparison({
                      refSampleName: e.target.value,
                    });
                  }}
                  custom
                >
                  {nameOptions.map((name, index) => {
                    return (
                      <option key={index} value={name}>
                        {name}
                      </option>
                    );
                  })}
                </Control>
              </Group>
            </Col>
            <Col sm="4">
              <Group controlId="signatureSet">
                <Label>Reference Signature Set</Label>
                <Control
                  disabled={!refSignatureSetOptions.length}
                  as="select"
                  value={refSignatureSet}
                  onChange={(e) => {
                    dispatchProfileComparison({
                      refSignatureSet: e.target.value,
                    });
                  }}
                  custom
                >
                  <option value="0">Select</option>
                  {refSignatureSetOptions.map((signatureSet, index) => {
                    return (
                      <option key={index} value={signatureSet}>
                        {signatureSet}
                      </option>
                    );
                  })}
                </Control>
              </Group>
            </Col>
            <Col sm="2">
              <Group controlId="signatureSet">
                <Label>Compare</Label>
                <Control
                  value={refCompare}
                  onChange={(e) => {
                    dispatchProfileComparison({
                      refCompare: e.target.value,
                    });
                  }}
                ></Control>
              </Group>
            </Col>
            <Col sm="2" className="m-auto">
              <Button
                variant="primary"
                onClick={() =>
                  calculateR('profileComparisonRefSig', {
                    profileType: refProfileType,
                    sampleName: refSampleName,
                    signatureSet: refSignatureSet,
                    compare: refCompare,
                  })
                }
              >
                Calculate
              </Button>
            </Col>
          </Row>

          <div
            id="refPlotPath"
            style={{ display: refPlotURL.length ? 'block' : 'none' }}
          >
            <div className="d-flex">
              <a
                className="px-2 py-1"
                href={refPlotURL}
                download={refPlotURL.split('/').slice(-1)[0]}
              >
                Download Plot
              </a>
            </div>
            <div className="p-2 border rounded">
              <Row>
                <Col>
                  <img className="w-100 my-4" src={refPlotURL}></img>
                </Col>
              </Row>
            </div>
          </div>
        </div>
      </Form>

      <pre className="border rounded p-1 mt-5">
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
      </pre>
    </div>
  );
}
