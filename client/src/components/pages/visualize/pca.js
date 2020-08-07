import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchVisualizeResults,
  dispatchPCA,
} from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

const { Group, Label, Control } = Form;

export default function PCA({ downloadResults, submitR, getRefSigOptions }) {
  const { pcWithinURL, pcRefSigURL } = useSelector(
    (state) => state.visualizeResults
  );
  const { profileOptions } = useSelector((state) => state.pyTab);
  const rootURL = window.location.pathname;
  const {
    profileType,
    signatureSet,
    signatureSetOptions,
    eig,
    pca1,
    pca2,
    heatmap,
    pca1Data,
    pca2Data,
    heatmapData,
    eigURL,
    pca1URL,
    pca2URL,
    heatmapURL,
    displayPCA,
    debugR,
    displayDebug,
    submitOverlay,
  } = useSelector((state) => state.pca);

  // load r plots after they are recieved
  useEffect(() => {
    if (eig && eig.length) {
      setRPlot(eig, 'eig');
    }
  }, [eig]);

  useEffect(() => {
    if (pca1 && pca1.length) {
      setRPlot(pca1, 'pca1');
    }
  }, [pca1]);

  useEffect(() => {
    if (pca2 && pca2.length) {
      setRPlot(pca2, 'pca2');
    }
  }, [pca2]);

  useEffect(() => {
    if (heatmap && heatmap.length) {
      setRPlot(heatmap, 'heatmap');
    }
  }, [heatmap]);

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
  // }, [withinPlotPath, refPlotPath]);

  async function setRPlot(plotPath, type) {
    try {
      const response = await getRefSigOptions(profileType);

      if (!response.ok) {
        const { msg } = await response.json();
        dispatchError(msg);
      } else {
        const pic = await response.blob();
        const objectURL = URL.createObjectURL(pic);

        if (type == 'eig') {
          if (eigURL.length) URL.revokeObjectURL(eigURL);
          dispatchPCA({
            eigURL: objectURL,
          });
        } else if (type == 'pca1') {
          if (pca1URL.length) URL.revokeObjectURL(pca1URL);
          dispatchPCA({
            pca1URL: objectURL,
          });
        } else if (type == 'pca2') {
          if (pca2URL.length) URL.revokeObjectURL(pca2URL);
          dispatchPCA({
            pca2URL: objectURL,
          });
        } else if (type == 'heatmap') {
          if (heatmapURL.length) URL.revokeObjectURL(heatmapURL);
          dispatchPCA({
            heatmapURL: objectURL,
          });
        }
      }
    } catch (err) {
      dispatchError(err);
    }
  }

  // get Signature Reference Sets for dropdown options
  async function getSignatureSet(profileType) {
    if (profileType && profileType.length) {
      dispatchPCA({ submitOverlay: true });
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

          dispatchPCA({
            profileType: profileType,
            signatureSetOptions: signatureSetOptions,
            signatureSet: signatureSetOptions[0],
            submitOverlay: false,
          });
        } else {
          dispatchError(await response.json());
          dispatchPCA({ submitOverlay: false });
        }
      } catch (err) {
        dispatchError(err);
        dispatchPCA({ submitOverlay: false });
      }
    }
  }

  async function calculateR(fn, args) {
    dispatchPCA({ submitOverlay: true, debugR: '' });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();
        dispatchPCA({
          debugR: err,
          submitOverlay: false,
        });
      } else {
        const { debugR, output } = await response.json();

        dispatchPCA({
          eig: '',
          pca1: '',
          pca2: '',
          heatmap: '',
          pca1Data: '',
          pca2Data: '',
          heatmapData: '',
        });
        dispatchPCA({
          debugR: debugR,
          submitOverlay: false,
          eig: output.eig,
          pca1: output.pca1,
          pca2: output.pca2,
          heatmap: output.heatmap,
          pca1Data: output.pca1Data,
          pca2Data: output.pca2Data,
          heatmapData: output.heatmapData,
        });
      }
    } catch (err) {
      dispatchError(err);
      dispatchPCA({ submitOverlay: false });
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
              dispatchPCA({
                displayPCA: !displayPCA,
              })
            }
          >
            Principal Component Analysis
          </Button>
        </Label>
        <div
          className="border rounded p-2"
          style={{ display: displayPCA ? 'block' : 'none' }}
        >
          <Row className="justify-content-center">
            <Col sm="5">
              <Group controlId="profileType">
                <Label>Profile Type</Label>
                <Control
                  as="select"
                  value={profileType}
                  onChange={(e) => getSignatureSet(e.target.value)}
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
                  onChange={(e) => {
                    dispatchPCA({
                      signatureSet: e.target.value,
                    });
                  }}
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
                  calculateR('pca', {
                    profileType: profileType,
                    signatureSet: signatureSet,
                  })
                }
              >
                Calculate
              </Button>
            </Col>
          </Row>

          <div
            id="eigPlot"
            className="my-4"
            style={{ display: eigURL.length ? 'block' : 'none' }}
          >
            <div className="d-flex">
              <a
                className="px-2 py-1"
                href={eigURL}
                download={eigURL.split('/').slice(-1)[0]}
              >
                Download Plot
              </a>
            </div>
            <div className="p-2 border rounded">
              <Row>
                <Col>
                  <img className="w-100 my-4" src={eigURL}></img>
                </Col>
              </Row>
            </div>
          </div>

          <div
            id="pca1Plot"
            className="my-4"
            style={{ display: pca1URL.length ? 'block' : 'none' }}
          >
            <div className="d-flex">
              <a
                className="px-2 py-1"
                href={pca1URL}
                download={pca1URL.split('/').slice(-1)[0]}
              >
                Download Plot
              </a>
              <span className="ml-auto">
                <Button
                  className="px-2 py-1"
                  variant="link"
                  onClick={() => downloadResults(pca1Data)}
                >
                  Download Results
                </Button>
              </span>
            </div>
            <div className="p-2 border rounded">
              <Row>
                <Col>
                  <img className="w-100 my-4" src={pca1URL}></img>
                </Col>
              </Row>
            </div>
          </div>

          <div
            id="withinPlotPath"
            className="my-4"
            style={{ display: pca2URL.length ? 'block' : 'none' }}
          >
            <div className="d-flex">
              <a
                className="px-2 py-1"
                href={pca2URL}
                download={pca2URL.split('/').slice(-1)[0]}
              >
                Download Plot
              </a>
              <span className="ml-auto">
                <Button
                  className="px-2 py-1"
                  variant="link"
                  onClick={() => downloadResults(pca2Data)}
                >
                  Download Results
                </Button>
              </span>
            </div>
            <div className="p-2 border rounded">
              <Row>
                <Col>
                  <img className="w-100 my-4" src={pca2URL}></img>
                </Col>
              </Row>
            </div>
          </div>

          <div
            id="withinPlotPath"
            className="my-4"
            style={{ display: heatmapURL.length ? 'block' : 'none' }}
          >
            <div className="d-flex">
              <a
                className="px-2 py-1"
                href={heatmapURL}
                download={heatmapURL.split('/').slice(-1)[0]}
              >
                Download Plot
              </a>
              <span className="ml-auto">
                <Button
                  className="px-2 py-1"
                  variant="link"
                  onClick={() => downloadResults(heatmapData)}
                >
                  Download Results
                </Button>
              </span>
            </div>
            <div className="p-2 border rounded">
              <Row>
                <Col>
                  <img className="w-100 my-4" src={heatmapURL}></img>
                </Col>
              </Row>
            </div>
          </div>
        </div>
      </Form>

      <Button
        variant="link"
        className="p-0 mt-5"
        onClick={() =>
          dispatchPCA({
            displayDebug: !displayDebug,
          })
        }
      >
        R Debug
      </Button>
      <pre
        className="border rounded p-1 "
        style={{ display: displayDebug ? 'block' : 'none' }}
      >
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
