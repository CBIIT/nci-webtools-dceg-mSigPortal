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
  const { pcWithinURL, pcRefSigURL } = useSelector(
    (state) => state.visualizeResults
  );
  const { nameOptions, profileOptions } = useSelector((state) => state.pyTab);
  const rootURL = window.location.pathname;
  const {
    within,
    refSig,
    withinPlotPath,
    refSigPlotPath,
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
    if (refSigPlotPath && refSigPlotPath.length) {
      setRPlot(refSigPlotPath, 'refsig');
    }
  }, [refSigPlotPath]);

  // call r wrapper on load
  // useEffect(() => {
  //   if (!pcWithinURL.length && !pcRefSigURL.length) {
  //     calculateR('cosineSimilarityWithin', {
  //       profileType: profileType1,
  //       matrixSize: matrixSize.replace('-', ''),
  //     });
  //     calculateR('cosineSimilarityRefSig', {
  //       profileType: profileType2,
  //       signatureSet: signatureSet,
  //     });
  //   }
  // }, [withinPlotPath, refSigPlotPath]);

  // retrieve signature set options on change
  useEffect(() => {
    // getSignatureSet(refSig.profileType);
  }, [refSig.profileType]);

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
          if (pcWithinURL.length) URL.revokeObjectURL(pcWithinURL);
          dispatchVisualizeResults({
            pcWithinURL: objectURL,
          });
        } else {
          if (pcRefSigURL.length) URL.revokeObjectURL(pcRefSigURL);
          dispatchVisualizeResults({
            pcRefSigURL: objectURL,
          });
        }
      }
    } catch (err) {
      dispatchError(err);
    }
  }

  // get Signature Reference Sets for dropdown options
  async function getSignatureSet(profileType) {
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
          refSig: {
            ...refSig,
            signatureSetOptions: signatureSetOptions,
            signatureSet: signatureSetOptions[0],
          },
          submitOverlay: false,
        });
        console.log('done');
      } else {
        dispatchError(await response.json());
        dispatchProfileComparison({ submitOverlay: false });
      }
    } catch (err) {
      dispatchError(err);
      dispatchProfileComparison({ submitOverlay: false });
    }
  }

  async function calculateR(fn, args) {
    dispatchProfileComparison({ submitOverlay: true });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();
        dispatchProfileComparison({
          debugR: err,
          submitOverlay: false,
        });
      } else {
        const data = await response.json();
        let update = {
          debugR: data.debugR,
          submitOverlay: false,
        };
        if (fn == 'profileComparisonWithin')
          update = {
            ...update,
            withinPlotPath: data.plot,
          };
        else {
          update = {
            ...update,
            refSigPlotPath: data.plot,
          };
        }
        dispatchProfileComparison(update);
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
                  value={within.profileType}
                  onChange={(e) =>
                    dispatchProfileComparison({
                      within: { ...within, profileType: e.target.value },
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
                value={within.sampleName1}
                onChange={(e) =>
                  dispatchProfileComparison({
                    within: { ...within, sampleName1: e.target.value },
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
                value={within.sampleName2}
                onChange={(e) =>
                  dispatchProfileComparison({
                    within: { ...within, sampleName2: e.target.value },
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
                    profileType: within.profileType,
                    sampleName1: within.sampleName1,
                    sampleName2: within.sampleName2,
                  })
                }
              >
                Calculate
              </Button>
            </Col>
          </Row>

          <div
            id="pcWithinPlot"
            style={{ display: pcWithinURL.length ? 'block' : 'none' }}
          >
            <div className="d-flex">
              <a
                className="px-2 py-1"
                href={pcWithinURL}
                download={pcWithinURL.split('/').slice(-1)[0]}
              >
                Download Plot
              </a>
            </div>
            <div className="mt-2 p-2 border rounded">
              <Row>
                <Col>
                  <img className="w-100 my-4" src={pcWithinURL}></img>
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
                  value={refSig.profileType}
                  onChange={(e) => {
                    dispatchProfileComparison({
                      refSig: { ...refSig, profileType: e.target.value },
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
                  value={refSig.sampleName}
                  onChange={(e) => {
                    dispatchProfileComparison({
                      refSig: { ...refSig, sampleName: e.target.value },
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
                  disabled={!refSig.signatureSetOptions.length}
                  as="select"
                  value={refSig.signatureSet}
                  onChange={(e) => {
                    dispatchProfileComparison({
                      refSig: { ...refSig, signatureSet: e.target.value },
                    });
                  }}
                  custom
                >
                  <option value="0">Select</option>
                  {refSig.signatureSetOptions.map((signatureSet, index) => {
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
                  value={refSig.compare}
                  onChange={(e) => {
                    dispatchProfileComparison({
                      refSig: { ...refSig, compare: e.target.value },
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
                    profileType: refSig.profileType,
                    sampleName: refSig.sampleName,
                    signatureSet: refSig.signatureSet,
                    compare: refSig.compare,
                  })
                }
              >
                Calculate
              </Button>
            </Col>
          </Row>

          <div
            id="refSigPlotPath"
            style={{ display: pcRefSigURL.length ? 'block' : 'none' }}
          >
            <div className="d-flex">
              <a
                className="px-2 py-1"
                href={pcRefSigURL}
                download={pcRefSigURL.split('/').slice(-1)[0]}
              >
                Download Plot
              </a>
            </div>
            <div className="mt-2 p-2 border rounded">
              <Row>
                <Col>
                  <img className="w-100 my-4" src={pcRefSigURL}></img>
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
