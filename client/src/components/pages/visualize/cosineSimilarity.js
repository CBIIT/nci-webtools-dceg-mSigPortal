import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchVisualizeResults,
  dispatchCosineSimilarity,
} from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

const { Group, Label, Control } = Form;

export default function CosineSimilarity({ submitR }) {
  const { csWithinURL, csRefSigURL } = useSelector(
    (state) => state.visualizeResults
  );
  const { mapping, profileOptions } = useSelector((state) => state.pyTab);
  const rootURL = window.location.pathname;
  const {
    profileType1,
    matrixSize,
    matrixOptions,
    profileType2,
    signatureSet,
    signatureSetOptions,
    csWithinPlot,
    csWithinTxt,
    csRefSigPlot,
    csRefSigTxt,
    submitOverlay,
    displayWithin,
    displayRefSig,
    debugR,
  } = useSelector((state) => state.cosineSimilarity);

  // load r plots after they are recieved
  useEffect(() => {
    if (csWithinPlot && csWithinPlot.length) {
      setRPlot(csWithinPlot, 'within');
    }
  }, [csWithinPlot]);

  useEffect(() => {
    if (csRefSigPlot && csRefSigPlot.length) {
      setRPlot(csRefSigPlot, 'refsig');
    }
  }, [csRefSigPlot]);

  // call r wrapper on load
  // useEffect(() => {
  //   if (!csWithinURL.length && !csRefSigURL.length) {
  //     calculateR('cosineSimilarityWithin', {
  //       profileType: profileType1,
  //       matrixSize: matrixSize.replace('-', ''),
  //     });
  //     calculateR('cosineSimilarityRefSig', {
  //       profileType: profileType2,
  //       signatureSet: signatureSet,
  //     });
  //   }
  // }, [csWithinPlot, csRefSigPlot]);

  // retrieve signature set options on change
  useEffect(() => {
    getSignatureSet(profileType2);
  }, [profileType2]);

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
          if (csWithinURL.length) URL.revokeObjectURL(csWithinURL);
          dispatchVisualizeResults({
            csWithinURL: objectURL,
          });
        } else {
          if (csRefSigURL.length) URL.revokeObjectURL(csRefSigURL);
          dispatchVisualizeResults({
            csRefSigURL: objectURL,
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
      dispatchCosineSimilarity({ submitOverlay: true });
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

          dispatchCosineSimilarity({
            signatureSetOptions: signatureSetOptions,
            signatureSet: signatureSetOptions[0],
            submitOverlay: false,
          });
        } else {
          dispatchError(await response.json());
          dispatchCosineSimilarity({ submitOverlay: false });
        }
      } catch (err) {
        dispatchError(err);
        dispatchCosineSimilarity({ submitOverlay: false });
      }
    }
  }

  async function calculateR(fn, args) {
    dispatchCosineSimilarity({ submitOverlay: true });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();
        dispatchCosineSimilarity({
          debugR: err,
          submitOverlay: false,
        });
      } else {
        const data = await response.json();
        let update = {
          debugR: data.debugR,
          submitOverlay: false,
        };
        if (fn == 'cosineSimilarityWithin')
          update = {
            ...update,
            csWithinPlot: data.plot,
            csWithinTxt: data.txt,
          };
        else {
          update = {
            ...update,
            csRefSigPlot: data.plot,
            csRefSigTxt: data.txt,
          };
        }
        dispatchCosineSimilarity(update);
      }
    } catch (err) {
      dispatchError(err);
      dispatchCosineSimilarity({ submitOverlay: false });
    }
  }

  //   download text results files
  async function downloadResults(txtPath) {
    try {
      const response = await fetch(`${rootURL}visualize/txt`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ path: txtPath }),
      });
      if (!response.ok) {
        const { msg } = await response.json();
        dispatchError(msg);
      } else {
        const file = await response.blob();
        const objectURL = URL.createObjectURL(file);
        const tempLink = document.createElement('a');

        tempLink.href = objectURL;
        tempLink.download = txtPath.split('/').slice(-1)[0];
        document.body.appendChild(tempLink);
        tempLink.click();
        document.body.removeChild(tempLink);
        URL.revokeObjectURL(objectURL);
      }
    } catch (err) {
      dispatchError(err);
    }
  }

  function handleProfileType1(profileType) {
    const matrixOptions = [
      ...new Set(
        mapping
          .filter((plot) => plot.Profile_Type == profileType)
          .map((plot) => plot.Matrix)
      ),
    ];

    dispatchCosineSimilarity({
      profileType1: profileType,
      matrixSize: matrixOptions[0],
      matrixOptions: matrixOptions,
    });
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
              dispatchCosineSimilarity({
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
            <Col sm="5">
              <Group controlId="profileType1">
                <Label>Profile Type</Label>
                <Control
                  as="select"
                  value={profileType1}
                  onChange={(e) => handleProfileType1(e.target.value)}
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
            <Col sm="5">
              <Label>Matrix Size</Label>
              <Control
                as="select"
                value={matrixSize}
                onChange={(e) =>
                  dispatchCosineSimilarity({
                    matrixSize: e.target.value,
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
            <Col sm="2" className="m-auto">
              <Button
                variant="primary"
                onClick={() =>
                  calculateR('cosineSimilarityWithin', {
                    profileType: profileType1,
                    matrixSize: matrixSize.replace('-', ''),
                  })
                }
              >
                Calculate
              </Button>
            </Col>
          </Row>

          <div
            id="csWithinPlot"
            style={{ display: csWithinURL.length ? 'block' : 'none' }}
          >
            <div className="d-flex">
              <a
                className="px-2 py-1"
                href={csWithinURL}
                download={csWithinURL.split('/').slice(-1)[0]}
              >
                Download Plot
              </a>
              <span className="ml-auto">
                <Button
                  className="px-2 py-1"
                  variant="link"
                  onClick={() => downloadResults(csWithinTxt)}
                >
                  Download Results
                </Button>
              </span>
            </div>
            <div className="mt-2 p-2 border rounded">
              <Row>
                <Col>
                  <img className="w-100 my-4" src={csWithinURL}></img>
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
              dispatchCosineSimilarity({
                displayRefSig: !displayRefSig,
              })
            }
          >
            Cosine Similarity to Reference Signatures
          </Button>
        </Label>
        <div
          className="border rounded p-2"
          style={{ display: displayRefSig ? 'block' : 'none' }}
        >
          <Row className="justify-content-center">
            <Col sm="5">
              <Group controlId="profileType2">
                <Label>Profile Type</Label>
                <Control
                  as="select"
                  value={profileType2}
                  onChange={(e) => {
                    dispatchCosineSimilarity({
                      profileType2: e.target.value,
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
            <Col sm="5">
              <Group controlId="signatureSet">
                <Label>Reference Signature Set</Label>
                <Control
                  disabled={!signatureSetOptions.length}
                  as="select"
                  value={signatureSet}
                  onChange={(e) =>
                    dispatchCosineSimilarity({
                      signatureSet: e.target.value,
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
            <Col sm="2" className="m-auto">
              <Button
                variant="primary"
                onClick={() =>
                  calculateR('cosineSimilarityRefSig', {
                    profileType: profileType2,
                    signatureSet: signatureSet,
                  })
                }
              >
                Calculate
              </Button>
            </Col>
          </Row>

          <div
            id="csRefSigPlot"
            style={{ display: csRefSigURL.length ? 'block' : 'none' }}
          >
            <div className="d-flex">
              <a
                className="px-2 py-1"
                href={csRefSigURL}
                download={csRefSigURL.split('/').slice(-1)[0]}
              >
                Download Plot
              </a>
              <span className="ml-auto">
                <Button
                  className="px-2 py-1"
                  variant="link"
                  onClick={() => downloadResults(csRefSigTxt)}
                >
                  Download Results
                </Button>
              </span>
            </div>
            <div className="mt-2 p-2 border rounded">
              <Row>
                <Col>
                  <img className="w-100 my-4" src={csRefSigURL}></img>
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
