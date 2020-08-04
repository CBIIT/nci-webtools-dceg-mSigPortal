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

export default function CosineSimilarity({ downloadResults, submitR }) {
  const { csWithinURL, csRefSigURL } = useSelector(
    (state) => state.visualizeResults
  );
  const { mapping, profileOptions } = useSelector((state) => state.pyTab);
  const rootURL = window.location.pathname;
  const {
    withinProfileType,
    withinMatrixSize,
    withinMatrixOptions,
    refProfileType,
    refSignatureSet,
    refSignatureSetOptions,
    withinPlotPath,
    withinTxtPath,
    refPlotPath,
    refTxtPath,
    submitOverlay,
    displayWithin,
    displayRefSig,
    debugR,
  } = useSelector((state) => state.cosineSimilarity);

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
  //   if (!csWithinURL.length && !csRefSigURL.length) {
  //     calculateR('cosineSimilarityWithin', {
  //       profileType: withinProfileType,
  //       matrixSize: withinMatrixSize.replace('-', ''),
  //     });
  //     calculateR('cosineSimilarityRefSig', {
  //       profileType: refProfileType,
  //       signatureSet: refSignatureSet,
  //     });
  //   }
  // }, [withinPlotPath, refPlotPath]);

  // retrieve signature set options on change
  useEffect(() => {
    getrefSignatureSet(refProfileType);
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
  async function getrefSignatureSet(profileType) {
    if (profileType && profileType.length) {
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
          const refSignatureSetOptions = await response.json();

          dispatchCosineSimilarity({
            refSignatureSetOptions: refSignatureSetOptions,
            refSignatureSet: refSignatureSetOptions[0],
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
        const { debugR, output } = await response.json();
        let update = {
          debugR: debugR,
          submitOverlay: false,
        };
        if (fn == 'cosineSimilarityWithin')
          dispatchCosineSimilarity({
            ...update,
            withinPlotPath: output.plotPath,
            withinTxtPath: output.txtPath,
          });
        else {
          dispatchCosineSimilarity({
            ...update,
            refPlotPath: output.plotPath,
            refTxtPath: output.txtPath,
          });
        }
      }
    } catch (err) {
      dispatchError(err);
      dispatchCosineSimilarity({ submitOverlay: false });
    }
  }

 

  function handlewithinProfileType(profileType) {
    const withinMatrixOptions = [
      ...new Set(
        mapping
          .filter((plot) => plot.Profile_Type == profileType)
          .map((plot) => plot.Matrix)
      ),
    ];

    dispatchCosineSimilarity({
      withinProfileType: profileType,
      withinMatrixSize: withinMatrixOptions[0],
      withinMatrixOptions: withinMatrixOptions,
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
              <Group controlId="withinProfileType">
                <Label>Profile Type</Label>
                <Control
                  as="select"
                  value={withinProfileType}
                  onChange={(e) => handlewithinProfileType(e.target.value)}
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
                value={withinMatrixSize}
                onChange={(e) =>
                  dispatchCosineSimilarity({
                    withinMatrixSize: e.target.value,
                  })
                }
                custom
              >
                <option value="0" disabled>
                  Select
                </option>
                {withinMatrixOptions.map((matrix, index) => {
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
                    profileType: withinProfileType,
                    matrixSize: withinMatrixSize.replace('-', ''),
                  })
                }
              >
                Calculate
              </Button>
            </Col>
          </Row>

          <div
            id="withinPlotPath"
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
                  onClick={() => downloadResults(withinTxtPath)}
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
              <Group controlId="refProfileType">
                <Label>Profile Type</Label>
                <Control
                  as="select"
                  value={refProfileType}
                  onChange={(e) => {
                    dispatchCosineSimilarity({
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
            <Col sm="5">
              <Group controlId="refSignatureSet">
                <Label>Reference Signature Set</Label>
                <Control
                  disabled={!refSignatureSetOptions.length}
                  as="select"
                  value={refSignatureSet}
                  onChange={(e) =>
                    dispatchCosineSimilarity({
                      refSignatureSet: e.target.value,
                    })
                  }
                  custom
                >
                  <option value="0">Select</option>
                  {refSignatureSetOptions.map((refSignatureSet, index) => {
                    return (
                      <option key={index} value={refSignatureSet}>
                        {refSignatureSet}
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
                    profileType: refProfileType,
                    signatureSet: refSignatureSet,
                  })
                }
              >
                Calculate
              </Button>
            </Col>
          </Row>

          <div
            id="refPlotPath"
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
                  onClick={() => downloadResults(refTxtPath)}
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
