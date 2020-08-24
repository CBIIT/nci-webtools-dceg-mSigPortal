import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { dispatchError, dispatchPCA } from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

const { Group, Label, Control } = Form;

export default function PCA({ downloadResults, submitR, getRefSigOptions }) {
  const { displayTab } = useSelector((state) => state.visualizeResults);
  const { profileOptions } = useSelector((state) => state.mutationalProfiles);
  const rootURL = window.location.pathname;
  const {
    profileType,
    signatureSet,
    signatureSetOptions,
    pca1,
    pca2,
    pca3,
    heatmap,
    pca2Data,
    pca3Data,
    heatmapData,
    pca1URL,
    pca2URL,
    pca3URL,
    heatmapURL,
    displayPCA,
    debugR,
    displayDebug,
    submitOverlay,
  } = useSelector((state) => state.pca);

  // load r plots after they are recieved
  // useEffect(() => {
  //   if (displayTab == 'pca') {
  //     if (pca1 && !pca1URL) setRPlot(pca1, 'pca1');
  //     if (pca2 && !pca2URL) setRPlot(pca2, 'pca2');
  //     if (pca3 && !pca3URL) setRPlot(pca3, 'pca3');
  //     if (heatmap && !heatmapURL) setRPlot(heatmap, 'heatmap');
  //   }
  // }, [pca1, pca2, pca3, heatmap, displayTab]);

  // useEffect(() => {
  //   if (
  //     profileType.length &&
  //     signatureSet.length &&
  //     !pca1 &&
  //     !pca2 &&
  //     !pca3 &&
  //     !heatmap &&
  //     !submitOverlay &&
  //     displayTab == 'pca'
  //   ) {
  //     calculateR('pca', {
  //       profileType: profileType,
  //       signatureSet: signatureSet,
  //     });
  //   }
  // }, [displayTab]);

  async function setRPlot(plotPath, type) {
    try {
      const response = await fetch(`${rootURL}visualize/svg`, {
        method: 'POST',
        headers: {
          Accept: 'image/svg',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ path: plotPath }),
      });

      if (!response.ok) {
        // console.log(await response.json());
      } else {
        const pic = await response.blob();
        const objectURL = URL.createObjectURL(pic);

        if (type == 'pca1') {
          if (pca1URL.length) URL.revokeObjectURL(pca1URL);
          dispatchPCA({
            pca1URL: objectURL,
          });
        } else if (type == 'pca2') {
          if (pca2URL.length) URL.revokeObjectURL(pca2URL);
          dispatchPCA({
            pca2URL: objectURL,
          });
        } else if (type == 'pca3') {
          if (pca3URL.length) URL.revokeObjectURL(pca3URL);
          dispatchPCA({
            pca3URL: objectURL,
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
        const response = await getRefSigOptions(profileType);

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
          debugR: debugR,
          submitOverlay: false,
          pca1: output.pca1,
          pca2: output.pca2,
          pca3: output.pca3,
          heatmap: output.heatmap,
          pca2Data: output.pca2Data,
          pca3Data: output.pca3Data,
          heatmapData: output.heatmapData,
        });
        setRPlot(output.pca1, 'pca1');
        setRPlot(output.pca2, 'pca2');
        setRPlot(output.pca3, 'pca3');
        setRPlot(output.heatmap, 'heatmap');
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
            </div>
            <div className="p-2 border rounded">
              <Row>
                <Col>
                  <img className="w-100 my-4 h-500" src={pca1URL}></img>
                </Col>
              </Row>
            </div>
          </div>

          <div
            id="pca2Plot"
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
                  <img className="w-100 my-4 h-600" src={pca2URL}></img>
                </Col>
              </Row>
            </div>
          </div>

          <div
            id="withinPlotPath"
            className="my-4"
            style={{ display: pca3URL.length ? 'block' : 'none' }}
          >
            <div className="d-flex">
              <a
                className="px-2 py-1"
                href={pca3URL}
                download={pca3URL.split('/').slice(-1)[0]}
              >
                Download Plot
              </a>
              <span className="ml-auto">
                <Button
                  className="px-2 py-1"
                  variant="link"
                  onClick={() => downloadResults(pca3Data)}
                >
                  Download Results
                </Button>
              </span>
            </div>
            <div className="p-2 border rounded">
              <Row>
                <Col>
                  <img className="w-100 my-4 h-600" src={pca3URL}></img>
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
                  <img className="w-100 my-4 h-600" src={heatmapURL}></img>
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
