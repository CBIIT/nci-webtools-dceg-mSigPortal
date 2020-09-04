import React, { useEffect } from 'react';
import { Form, Row, Col, Button, Accordion, Card } from 'react-bootstrap';
import Select from 'react-select';
import { useSelector } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import { dispatchError, dispatchPCA } from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

const { Group, Label } = Form;
const { Header, Body } = Card;
const { Toggle, Collapse } = Accordion;

export default function PCA({ downloadResults, submitR, getRefSigOptions }) {
  const { matrixList } = useSelector((state) => state.visualizeResults);
  const {
    source,
    study,
    studyOptions,
    cancerType,
    pExperimentalStrategy,
    pDataOptions,
  } = useSelector((state) => state.visualize);
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
    pca1Err,
    pca2Err,
    pca3Err,
    heatmapErr,
    debugR,
    displayDebug,
    submitOverlay,
    pubPCA,

    userProfileType,
    userMatrixSize,
    userMatrixOptions,
    pubStudy,
    pubCancerType,
    pubCancerTypeOptions,
    pubPca1,
    pubPca2,
    pubPca3,
    pubPca2Data,
    pubPca3Data,
    pubPca1URL,
    pubPca2URL,
    pubPca3URL,
    displayPub,
    pubPca1Err,
    pubPca2Err,
    pubPca3Err,
    pubSubmitOverlay,
  } = useSelector((state) => state.pca);

  const selectFix = {
    styles: {
      menuPortal: (base) => ({ ...base, zIndex: 9999 }),
    },
    menuPortalTarget: document.body,
  };

  async function setRPlot(plotPath, type) {
    if (plotPath) {
      dispatchPCA({ submitOverlay: true });

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
            if (pca1URL) URL.revokeObjectURL(pca1URL);
            dispatchPCA({ pca1URL: objectURL });
          } else if (type == 'pca2') {
            if (pca2URL) URL.revokeObjectURL(pca2URL);
            dispatchPCA({ pca2URL: objectURL });
          } else if (type == 'pca3') {
            if (pca3URL) URL.revokeObjectURL(pca3URL);
            dispatchPCA({ pca3URL: objectURL });
          } else if (type == 'heatmap') {
            if (heatmapURL) URL.revokeObjectURL(heatmapURL);
            dispatchPCA({ heatmapURL: objectURL });
          } else if (type == 'pubPca1') {
            if (pubPca1URL) URL.revokeObjectURL(pubPca1);
            dispatchPCA({ pubPca1URL: objectURL });
          } else if (type == 'pubPca2') {
            if (pubPca2URL) URL.revokeObjectURL(pubPca2URL);
            dispatchPCA({ pubPca2URL: objectURL });
          } else if (type == 'pubPca3') {
            if (pubPca3URL) URL.revokeObjectURL(pubPca3URL);
            dispatchPCA({ pubPca3URL: objectURL });
          }
        }
      } catch (err) {
        dispatchError(err);
      }
      dispatchPCA({ submitOverlay: false });
    } else {
      if (type == 'pca1') {
        if (pca1URL) URL.revokeObjectURL(pca1URL);
        dispatchPCA({ pca1Err: true, pca1URL: '' });
      } else if (type == 'pca2') {
        if (pca2URL) URL.revokeObjectURL(pca2URL);
        dispatchPCA({ pca2Err: true, pca2URL: '' });
      } else if (type == 'pca3') {
        if (pca3URL) URL.revokeObjectURL(pca3URL);
        dispatchPCA({ pca3Err: true, pca3URL: '' });
      } else if (type == 'heatmap') {
        if (heatmapURL) URL.revokeObjectURL(heatmapURL);
        dispatchPCA({ heatmapErr: true, heatmapURL: '' });
      } else if (type == 'pubPca1') {
        if (pubPca1URL) URL.revokeObjectURL(pubPca1URL);
        dispatchPCA({ pubPca1Err: true, pubPca1URL: '' });
      } else if (type == 'pubPca2') {
        if (pubPca2URL) URL.revokeObjectURL(pubPca2URL);
        dispatchPCA({ pubPca2Err: true, pubPca2URL: '' });
      } else if (type == 'pubPca3') {
        if (pubPca3URL) URL.revokeObjectURL(pubPca3URL);
        dispatchPCA({ pubPca3Err: true, pubPca3URL: '' });
      }
    }
  }

  // get Signature Reference Sets for dropdown options
  async function getSignatureSet(profileType) {
    if (profileType) {
      dispatchPCA({ submitOverlay: true });
      try {
        const response = await getRefSigOptions(profileType);

        if (response.ok) {
          const signatureSetOptions = await response.json();

          dispatchPCA({
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
    if (!fn.includes('With')) {
      dispatchPCA({
        submitOverlay: true,
        debugR: '',
        pca1Err: false,
        pca2Err: false,
        pca3Err: false,
        heatmapErr: false,
      });
    } else {
      dispatchPCA({
        pubSubmitOverlay: true,
        debugR: '',
        pubPca1Err: false,
        pubPca2Err: false,
        pubPca3Err: false,
      });
    }

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();

        dispatchPCA({ debugR: err });
        if (!fn.includes('With')) {
          dispatchPCA({ submitOverlay: false });
        } else {
          dispatchPCA({
            pubSubmitOverlay: false,
          });
        }
      } else {
        const { debugR, output } = await response.json();

        if (!fn.includes('With')) {
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
        } else {
          dispatchPCA({
            debugR: debugR,
            pubSubmitOverlay: false,
            pubPca1: output.pca1,
            pubPca2: output.pca2,
            pubPca3: output.pca3,

            pubPca2Data: output.pca2Data,
            pubPca3Data: output.pca3Data,
          });
          setRPlot(output.pca1, 'pubPca1');
          setRPlot(output.pca2, 'pubPca2');
          setRPlot(output.pca3, 'pubPca3');
        }
      }
    } catch (err) {
      dispatchError(err);
      if (!fn.includes('With')) {
        dispatchPCA({ submitOverlay: false });
      } else {
        dispatchPCA({
          pubSubmitOverlay: false,
        });
      }
    }
  }

  function handleProfileType(profileType) {
    const userMatrixOptions = [
      ...new Set(
        matrixList
          .filter((matrix) => matrix.Profile_Type == profileType)
          .map((matrix) => matrix.Matrix_Size)
      ),
    ];

    dispatchPCA({
      userProfileType: profileType,
      userMatrixSize: userMatrixOptions[0],
      userMatrixOptions: userMatrixOptions,
    });
  }

  function handleStudyChange(study) {
    const cancerTypeOptions = [
      ...new Set(
        pDataOptions
          .filter((data) => data.Study == study)
          .map((data) => data.Cancer_Type)
      ),
    ];

    dispatchPCA({
      pubStudy: study,
      pubCancerType: cancerTypeOptions[0],
      pubCancerTypeOptions: cancerTypeOptions,
    });
  }

  function handleCancerChange(cancer) {
    dispatchPCA({
      pubCancerType: cancer,
    });
  }

  return (
    <div>
      <Accordion defaultActiveKey="0">
        <Card>
          <Toggle
            className="font-weight-bold"
            as={Header}
            eventKey="0"
            onClick={() =>
              dispatchPCA({
                displayPCA: !displayPCA,
              })
            }
          >
            {displayPCA == true ? (
              <FontAwesomeIcon icon={faMinus} />
            ) : (
              <FontAwesomeIcon icon={faPlus} />
            )}{' '}
            Principal Component Analysis
          </Toggle>
          <Collapse eventKey="0">
            <Body>
              <Form>
                <LoadingOverlay active={submitOverlay} />
                <div>
                  <Row className="justify-content-center">
                    <Col sm="5">
                      <Group controlId="profileType">
                        <Label>Profile Type</Label>
                        <Select
                          options={profileOptions}
                          value={[profileType]}
                          onChange={(profileType) => {
                            dispatchPCA({
                              profileType: profileType,
                            });
                            getSignatureSet(profileType);
                          }}
                          getOptionLabel={(option) => option}
                          getOptionValue={(option) => option}
                          {...selectFix}
                        />
                      </Group>
                    </Col>

                    <Col sm="5">
                      <Group controlId="signatureSet">
                        <Label>Reference Signature Set</Label>

                        <Select
                          options={signatureSetOptions}
                          value={[signatureSet]}
                          onChange={(signatureSet) => {
                            dispatchPCA({
                              signatureSet: signatureSet,
                            });
                          }}
                          getOptionLabel={(option) => option}
                          getOptionValue={(option) => option}
                          {...selectFix}
                        />
                      </Group>
                    </Col>
                    <Col sm="1" className="m-auto">
                      <Button
                        variant="primary"
                        onClick={() => {
                          if (source == 'user') {
                            calculateR('pca', {
                              profileType: profileType,
                              signatureSet: signatureSet,
                              matrixList: JSON.stringify(
                                matrixList.filter(
                                  (matrix) => matrix.Profile_Type == profileType
                                )
                              ),
                            });
                          } else {
                            calculateR('pcaPublic', {
                              profileType: profileType,
                              signatureSet: signatureSet,
                              study: study,
                              cancerType: cancerType,
                              experimentalStrategy: pExperimentalStrategy,
                            });
                          }
                        }}
                      >
                        Calculate
                      </Button>
                    </Col>
                  </Row>

                  <div id="pca1Plot">
                    <div style={{ display: pca1Err ? 'block' : 'none' }}>
                      <h4>PCA 1</h4>
                      <p>
                        An error has occured. Check the debug section for more
                        info.
                      </p>
                    </div>
                    <div
                      className="my-4"
                      style={{ display: pca1URL ? 'block' : 'none' }}
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
                            <img
                              className="w-100 my-4 h-500"
                              src={pca1URL}
                            ></img>
                          </Col>
                        </Row>
                      </div>
                    </div>
                  </div>

                  <div id="pca2Plot">
                    <div style={{ display: pca2Err ? 'block' : 'none' }}>
                      <h4>PCA 2</h4>
                      <p>
                        An error has occured. Check the debug section for more
                        info.
                      </p>
                    </div>
                    <div
                      className="my-4"
                      style={{ display: pca2URL ? 'block' : 'none' }}
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
                            <img
                              className="w-100 my-4 h-600"
                              src={pca2URL}
                            ></img>
                          </Col>
                        </Row>
                      </div>
                    </div>
                  </div>

                  <div id="pca3Plot">
                    <div style={{ display: pca3Err ? 'block' : 'none' }}>
                      <h4>PCA 3</h4>
                      <p>
                        An error has occured. Check the debug section for more
                        info.
                      </p>
                    </div>
                    <div
                      className="my-4"
                      style={{ display: pca3URL ? 'block' : 'none' }}
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
                            <img
                              className="w-100 my-4 h-600"
                              src={pca3URL}
                            ></img>
                          </Col>
                        </Row>
                      </div>
                    </div>
                  </div>

                  <div id="heatmapPlot">
                    <div style={{ display: heatmapErr ? 'block' : 'none' }}>
                      <h4>Heatmap</h4>
                      <p>
                        An error has occured. Check the debug section for more
                        info.
                      </p>
                    </div>
                    <div
                      className="my-4"
                      style={{ display: heatmapURL ? 'block' : 'none' }}
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
                            <img
                              className="w-100 my-4 h-600"
                              src={heatmapURL}
                            ></img>
                          </Col>
                        </Row>
                      </div>
                    </div>
                  </div>
                </div>
              </Form>
            </Body>
          </Collapse>
        </Card>
      </Accordion>

      {source == 'user' && (
        <Accordion defaultActiveKey="1">
          <Card>
            <Toggle
              className="font-weight-bold"
              as={Header}
              eventKey="1"
              onClick={() => dispatchPCA({ displayPub: !displayPub })}
            >
              {displayPub == true ? (
                <FontAwesomeIcon icon={faMinus} />
              ) : (
                <FontAwesomeIcon icon={faPlus} />
              )}{' '}
              Principal Component Analysis with Public Data
            </Toggle>
            <Collapse eventKey="1">
              <Body>
                <Form>
                  <LoadingOverlay active={pubSubmitOverlay} />
                  <div>
                    <Row className="justify-content-center">
                      <Col sm="2">
                        <Group controlId="profileType">
                          <Label>Profile Type</Label>
                          <Select
                            options={profileOptions}
                            value={[userProfileType]}
                            onChange={(profile) => handleProfileType(profile)}
                            getOptionLabel={(option) => option}
                            getOptionValue={(option) => option}
                            {...selectFix}
                          />
                        </Group>
                      </Col>

                      <Col sm="2">
                        <Group controlId="signatureSet">
                          <Label>Matrix Size</Label>
                          <Select
                            options={userMatrixOptions}
                            value={[userMatrixSize]}
                            onChange={(matrix) => {
                              dispatchPCA({ userMatrixSize: matrix });
                            }}
                            getOptionLabel={(option) => option}
                            getOptionValue={(option) => option}
                            {...selectFix}
                          />
                        </Group>
                      </Col>

                      <Col sm="2">
                        <Group controlId="pubStudy">
                          <Label>Study</Label>
                          <Select
                            options={studyOptions}
                            value={[pubStudy]}
                            onChange={(study) => handleStudyChange(study)}
                            getOptionLabel={(option) => option}
                            getOptionValue={(option) => option}
                            {...selectFix}
                          />
                        </Group>
                      </Col>
                      <Col sm="4">
                        <Group controlId="pubCancerType">
                          <Label>Cancer Type</Label>
                          <Select
                            options={pubCancerTypeOptions}
                            value={[pubCancerType]}
                            onChange={(cancerType) =>
                              handleCancerChange(cancerType)
                            }
                            getOptionLabel={(option) => option}
                            getOptionValue={(option) => option}
                            {...selectFix}
                          />
                        </Group>
                      </Col>
                      <Col sm="1" className="m-auto">
                        <Button
                          variant="primary"
                          onClick={() =>
                            calculateR('pcaWithPublic', {
                              matrixFile: matrixList.filter(
                                (path) =>
                                  path.Profile_Type == userProfileType &&
                                  path.Matrix_Size == userMatrixSize
                              )[0].Path,
                              study: pubStudy,
                              cancerType: pubCancerType,
                              profileName: userProfileType + userMatrixSize,
                            })
                          }
                        >
                          Calculate
                        </Button>
                      </Col>
                    </Row>

                    <div id="pubPca1Plot">
                      <div style={{ display: pubPca1Err ? 'block' : 'none' }}>
                        <h4>PCA 1</h4>
                        <p>
                          An error has occured. Check the debug section for more
                          info.
                        </p>
                      </div>
                      <div
                        className="my-4"
                        style={{ display: pubPca1URL ? 'block' : 'none' }}
                      >
                        <div className="d-flex">
                          <a
                            className="px-2 py-1"
                            href={pubPca1URL}
                            download={pubPca1URL.split('/').slice(-1)[0]}
                          >
                            Download Plot
                          </a>
                        </div>
                        <div className="p-2 border rounded">
                          <Row>
                            <Col>
                              <img
                                className="w-100 my-4 h-500"
                                src={pubPca1URL}
                              ></img>
                            </Col>
                          </Row>
                        </div>
                      </div>
                    </div>

                    <div id="pubPca2Plot">
                      <div style={{ display: pubPca2Err ? 'block' : 'none' }}>
                        <h4>PCA 2</h4>
                        <p>
                          An error has occured. Check the debug section for more
                          info.
                        </p>
                      </div>
                      <div
                        className="my-4"
                        style={{ display: pubPca2URL ? 'block' : 'none' }}
                      >
                        <div className="d-flex">
                          <a
                            className="px-2 py-1"
                            href={pubPca2URL}
                            download={pubPca2URL.split('/').slice(-1)[0]}
                          >
                            Download Plot
                          </a>
                          <span className="ml-auto">
                            <Button
                              className="px-2 py-1"
                              variant="link"
                              onClick={() => downloadResults(pubPca2Data)}
                            >
                              Download Results
                            </Button>
                          </span>
                        </div>
                        <div className="p-2 border rounded">
                          <Row>
                            <Col>
                              <img
                                className="w-100 my-4 h-600"
                                src={pubPca2URL}
                              ></img>
                            </Col>
                          </Row>
                        </div>
                      </div>
                    </div>

                    <div id="pubPca3Plot">
                      <div style={{ display: pubPca3Err ? 'block' : 'none' }}>
                        <h4>PCA 3</h4>
                        <p>
                          An error has occured. Check the debug section for more
                          info.
                        </p>
                      </div>
                      <div
                        className="my-4"
                        style={{ display: pubPca3URL ? 'block' : 'none' }}
                      >
                        <div className="d-flex">
                          <a
                            className="px-2 py-1"
                            href={pubPca3URL}
                            download={pubPca3URL.split('/').slice(-1)[0]}
                          >
                            Download Plot
                          </a>
                          <span className="ml-auto">
                            <Button
                              className="px-2 py-1"
                              variant="link"
                              onClick={() => downloadResults(pubPca3Data)}
                            >
                              Download Results
                            </Button>
                          </span>
                        </div>
                        <div className="p-2 border rounded">
                          <Row>
                            <Col>
                              <img
                                className="w-100 my-4 h-600"
                                src={pubPca3URL}
                              ></img>
                            </Col>
                          </Row>
                        </div>
                      </div>
                    </div>
                  </div>
                </Form>
              </Body>
            </Collapse>
          </Card>
        </Accordion>
      )}

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
