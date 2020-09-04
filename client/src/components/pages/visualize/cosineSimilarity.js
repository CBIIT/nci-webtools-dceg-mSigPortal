import React, { useEffect } from 'react';
import { Form, Row, Col, Button, Accordion, Card } from 'react-bootstrap';
import Select from 'react-select';
import { useSelector } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import {
  dispatchError,
  dispatchCosineSimilarity,
} from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

const { Group, Label, Control } = Form;
const { Header, Body } = Card;
const { Toggle, Collapse } = Accordion;

export default function CosineSimilarity({
  downloadResults,
  submitR,
  getRefSigOptions,
}) {
  const {
    source,
    study,
    studyOptions,
    cancerType,
    pExperimentalStrategy,
    pDataOptions,
  } = useSelector((state) => state.visualize);
  const { profileOptions } = useSelector((state) => state.mutationalProfiles);
  const { matrixList, svgList } = useSelector(
    (state) => state.visualizeResults
  );
  const rootURL = window.location.pathname;
  const {
    withinProfileType,
    withinMatrixSize,
    withinMatrixOptions,
    refProfileType,
    refSignatureSet,
    refSignatureSetOptions,
    userProfileType,
    userMatrixSize,
    userMatrixOptions,
    pubStudy,
    pubCancerType,
    pubCancerTypeOptions,
    withinPlotPath,
    withinTxtPath,
    refPlotPath,
    refTxtPath,
    pubPlotPath,
    pubTxtPath,
    withinPlotURL,
    refPlotURL,
    pubPlotURL,
    displayWithin,
    displayRefSig,
    displayPublic,
    withinErr,
    refErr,
    pubErr,
    withinSubmitOverlay,
    refSubmitOverlay,
    pubSubmitOverlay,
    debugR,
    displayDebug,
  } = useSelector((state) => state.cosineSimilarity);

  const selectFix = {
    styles: {
      menuPortal: (base) => ({ ...base, zIndex: 9999 }),
    },
    menuPortalTarget: document.body,
  };

  function setOverlay(type, display) {
    if (type == 'within') {
      dispatchCosineSimilarity({ withinSubmitOverlay: display });
    } else if (type == 'refsig') {
      dispatchCosineSimilarity({ refSubmitOverlay: display });
    } else {
      dispatchCosineSimilarity({ pubSubmitOverlay: display });
    }
  }

  async function setRPlot(plotPath, type) {
    setOverlay(type, true);
    if (plotPath) {
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

          if (type == 'within') {
            if (withinPlotURL) URL.revokeObjectURL(withinPlotURL);
            dispatchCosineSimilarity({
              withinPlotURL: objectURL,
            });
          } else if (type == 'refsig') {
            if (refPlotURL) URL.revokeObjectURL(refPlotURL);
            dispatchCosineSimilarity({
              refPlotURL: objectURL,
            });
          } else {
            if (pubPlotURL) URL.revokeObjectURL(pubPlotURL);
            dispatchCosineSimilarity({
              pubPlotURL: objectURL,
            });
          }
          setOverlay(type, false);
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (type == 'within') {
        if (withinPlotURL) URL.revokeObjectURL(withinPlotURL);
        dispatchCosineSimilarity({ withinErr: true, withinPlotURL: '' });
      } else if (type == 'refsig') {
        if (refPlotURL) URL.revokeObjectURL(refPlotURL);
        dispatchCosineSimilarity({ refErr: true, refPlotURL: '' });
      } else {
        if (pubPlotURL) URL.revokeObjectURL(pubPlotURL);
        dispatchCosineSimilarity({ pubErr: true, pubPlotURL: '' });
      }
    }
    setOverlay(type, false);
  }

  // get Signature Reference Sets for dropdown options
  async function getSignatureSet(profileType) {
    if (profileType) {
      dispatchCosineSimilarity({ refSubmitOverlay: true });
      try {
        const response = await getRefSigOptions(profileType);

        if (response.ok) {
          const refSignatureSetOptions = await response.json();

          dispatchCosineSimilarity({
            refSignatureSetOptions: refSignatureSetOptions,
            refSignatureSet: refSignatureSetOptions[0],
            refSubmitOverlay: false,
          });
        } else {
          dispatchError(await response.json());
          dispatchCosineSimilarity({ refSubmitOverlay: false });
        }
      } catch (err) {
        dispatchError(err);
        dispatchCosineSimilarity({ refSubmitOverlay: false });
      }
    }
  }

  async function calculateR(fn, args) {
    if (fn.includes('cosineSimilarityWithin')) {
      console.log(fn);
      dispatchCosineSimilarity({
        withinSubmitOverlay: true,
        withinErr: false,
        debugR: '',
      });
    } else if (fn.includes('cosineSimilarityRefSig')) {
      dispatchCosineSimilarity({
        refSubmitOverlay: true,
        refErr: false,
        debugR: '',
      });
    } else {
      dispatchCosineSimilarity({
        pubSubmitOverlay: true,
        pubErr: false,
        debugR: '',
      });
    }
    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();
        if (fn.includes('cosineSimilarityWithin')) {
          dispatchCosineSimilarity({
            withinSubmitOverlay: false,
            debugR: err,
          });
        } else if (fn.includes('cosineSimilarityRefSig')) {
          dispatchCosineSimilarity({
            refSubmitOverlay: false,
            debugR: err,
          });
        } else {
          dispatchCosineSimilarity({
            pubSubmitOverlay: false,
            debugR: err,
          });
        }
      } else {
        const { debugR, output } = await response.json();

        if (fn.includes('cosineSimilarityWithin')) {
          dispatchCosineSimilarity({
            debugR: debugR,
            withinSubmitOverlay: false,
            withinPlotPath: output.plotPath,
            withinTxtPath: output.txtPath,
          });
          setRPlot(output.plotPath, 'within');
        } else if (fn.includes('cosineSimilarityRefSig')) {
          dispatchCosineSimilarity({
            debugR: debugR,
            refSubmitOverlay: false,
            refPlotPath: output.plotPath,
            refTxtPath: output.txtPath,
          });
          setRPlot(output.plotPath, 'refsig');
        } else {
          dispatchCosineSimilarity({
            debugR: debugR,
            pubSubmitOverlay: false,
            pubPlotPath: output.plotPath,
            pubTxtPath: output.txtPath,
          });
          setRPlot(output.plotPath, 'pub');
        }
      }
    } catch (err) {
      dispatchError(err);
      if (fn.includes('cosineSimilarityWithin')) {
        dispatchCosineSimilarity({ withinSubmitOverlay: false });
      } else if (fn.includes('cosineSimilarityRefSig')) {
        dispatchCosineSimilarity({ refSubmitOverlay: false });
      } else {
        dispatchCosineSimilarity({ pubSubmitOverlay: false });
      }
    }
  }

  function handleWithinProfileType(profileType) {
    if (source == 'user') {
      const withinMatrixOptions = [
        ...new Set(
          matrixList
            .filter((matrix) => matrix.Profile_Type == profileType)
            .map((matrix) => matrix.Matrix_Size)
        ),
      ];

      dispatchCosineSimilarity({
        withinProfileType: profileType,
        withinMatrixSize: withinMatrixOptions[0],
        withinMatrixOptions: withinMatrixOptions,
      });
    } else {
      const withinMatrixOptions = [
        ...new Set(
          svgList
            .filter((plot) => plot.Profile.indexOf(profileType) > -1)
            .map((plot) => plot.Profile.match(/\d+/gi)[0])
        ),
      ];

      dispatchCosineSimilarity({
        withinProfileType: profileType,
        withinMatrixSize: withinMatrixOptions[0],
        withinMatrixOptions: withinMatrixOptions,
      });
    }
  }

  function handlePublicProfileType(profileType) {
    const userMatrixOptions = [
      ...new Set(
        matrixList
          .filter((matrix) => matrix.Profile_Type == profileType)
          .map((matrix) => matrix.Matrix_Size)
      ),
    ];

    dispatchCosineSimilarity({
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

    dispatchCosineSimilarity({
      pubStudy: study,
      pubCancerType: cancerTypeOptions[0],
      pubCancerTypeOptions: cancerTypeOptions,
    });
  }

  function handleCancerChange(cancer) {
    dispatchCosineSimilarity({
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
              dispatchCosineSimilarity({
                displayWithin: !displayWithin,
              })
            }
          >
            {displayWithin == true ? (
              <FontAwesomeIcon icon={faMinus} />
            ) : (
              <FontAwesomeIcon icon={faPlus} />
            )}{' '}
            Cosine Similarity Within Samples
          </Toggle>
          <Collapse eventKey="0">
            <Body>
              <Form>
                <LoadingOverlay active={withinSubmitOverlay} />
                <div>
                  <Row className="justify-content-center">
                    <Col sm="5">
                      <Group controlId="withinProfileType">
                        <Label>Profile Type</Label>
                        <Select
                          options={profileOptions}
                          value={[withinProfileType]}
                          onChange={(profile) =>
                            handleWithinProfileType(profile)
                          }
                          getOptionLabel={(option) => option}
                          getOptionValue={(option) => option}
                          {...selectFix}
                        />
                      </Group>
                    </Col>
                    <Col sm="5">
                      <Label>Matrix Size</Label>
                      <Select
                        options={withinMatrixOptions}
                        value={[withinMatrixSize]}
                        onChange={(matrix) =>
                          dispatchCosineSimilarity({
                            withinMatrixSize: matrix,
                          })
                        }
                        getOptionLabel={(option) => option}
                        getOptionValue={(option) => option}
                        {...selectFix}
                      />
                    </Col>
                    <Col sm="1" className="m-auto">
                      <Button
                        variant="primary"
                        onClick={() => {
                          if (source == 'user') {
                            calculateR('cosineSimilarityWithin', {
                              matrixFile: matrixList.filter(
                                (path) =>
                                  path.Profile_Type == withinProfileType &&
                                  path.Matrix_Size == withinMatrixSize
                              )[0].Path,
                            });
                          } else {
                            calculateR('cosineSimilarityWithinPublic', {
                              profileType: withinProfileType,
                              matrixSize: withinMatrixSize,
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

                  <div id="withinPlot">
                    <div style={{ display: withinErr ? 'block' : 'none' }}>
                      <p>
                        An error has occured. Check the debug section for more
                        info.
                      </p>
                    </div>
                    <div style={{ display: withinPlotURL ? 'block' : 'none' }}>
                      <div className="d-flex">
                        <a
                          className="px-2 py-1"
                          href={withinPlotURL}
                          download={withinPlotURL.split('/').slice(-1)[0]}
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
                      <div className="p-2 border rounded">
                        <Row>
                          <Col>
                            <img
                              className="w-100 my-4 h-500"
                              src={withinPlotURL}
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

      <Accordion defaultActiveKey="1">
        <Card>
          <Toggle
            className="font-weight-bold"
            as={Header}
            eventKey="1"
            onClick={() =>
              dispatchCosineSimilarity({
                displayRefSig: !displayRefSig,
              })
            }
          >
            {displayRefSig == true ? (
              <FontAwesomeIcon icon={faMinus} />
            ) : (
              <FontAwesomeIcon icon={faPlus} />
            )}{' '}
            Cosine Similarity to Reference Signatures
          </Toggle>
          <Collapse eventKey="1">
            <Body>
              <Form className="my-2">
                <LoadingOverlay active={refSubmitOverlay} />
                <div>
                  <Row className="justify-content-center">
                    <Col sm="5">
                      <Group controlId="refProfileType">
                        <Label>Profile Type</Label>
                        <Select
                          options={profileOptions}
                          value={[refProfileType]}
                          onChange={(refProfileType) => {
                            dispatchCosineSimilarity({
                              refProfileType: refProfileType,
                            });
                            getSignatureSet(refProfileType);
                          }}
                          getOptionLabel={(option) => option}
                          getOptionValue={(option) => option}
                          {...selectFix}
                        />
                      </Group>
                    </Col>
                    <Col sm="5">
                      <Group controlId="refSignatureSet">
                        <Label>Reference Signature Set</Label>
                        <Select
                          options={refSignatureSetOptions}
                          value={[refSignatureSet]}
                          onChange={(refSignatureSet) => {
                            dispatchCosineSimilarity({
                              refSignatureSet: refSignatureSet,
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
                            calculateR('cosineSimilarityRefSig', {
                              profileType: refProfileType,
                              signatureSet: refSignatureSet,
                              matrixList: JSON.stringify(
                                matrixList.filter(
                                  (matrix) =>
                                    matrix.Profile_Type == refProfileType
                                )
                              ),
                            });
                          } else {
                            calculateR('cosineSimilarityRefSigPublic', {
                              profileType: refProfileType,
                              signatureSet: refSignatureSet,
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

                  <div id="refPlot">
                    <div style={{ display: refErr ? 'block' : 'none' }}>
                      <p>
                        An error has occured. Check the debug section for more
                        info.
                      </p>
                    </div>
                    <div style={{ display: refPlotURL ? 'block' : 'none' }}>
                      <div className="d-flex">
                        <a
                          className="px-2 py-1"
                          href={refPlotURL}
                          download={refPlotURL.split('/').slice(-1)[0]}
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
                      <div className="p-2 border rounded">
                        <Row>
                          <Col>
                            {refErr == true && (
                              <p>
                                An error has occured. Check the debug section
                                for more info.
                              </p>
                            )}
                            <img
                              className="w-100 my-4 h-500"
                              src={refPlotURL}
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
        <Accordion defaultActiveKey="2">
          <Card>
            <Toggle
              className="font-weight-bold"
              as={Header}
              eventKey="2"
              onClick={() =>
                dispatchCosineSimilarity({
                  displayPublic: !displayPublic,
                })
              }
            >
              {displayPublic == true ? (
                <FontAwesomeIcon icon={faMinus} />
              ) : (
                <FontAwesomeIcon icon={faPlus} />
              )}{' '}
              Cosine Similarity to Public Data
            </Toggle>
            <Collapse eventKey="2">
              <Body>
                <Form>
                  <LoadingOverlay active={pubSubmitOverlay} />
                  <div>
                    <Row className="justify-content-center">
                      <Col sm="2">
                        <Group controlId="userProfileType">
                          <Label>Profile Type</Label>
                          <Select
                            options={profileOptions}
                            value={[userProfileType]}
                            onChange={(profile) =>
                              handlePublicProfileType(profile)
                            }
                            getOptionLabel={(option) => option}
                            getOptionValue={(option) => option}
                            {...selectFix}
                          />
                        </Group>
                      </Col>
                      <Col sm="2">
                        <Label>Matrix Size</Label>
                        <Select
                          options={userMatrixOptions}
                          value={[userMatrixSize]}
                          onChange={(matrix) =>
                            dispatchCosineSimilarity({
                              userMatrixSize: matrix,
                            })
                          }
                          getOptionLabel={(option) => option}
                          getOptionValue={(option) => option}
                          {...selectFix}
                        />
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
                            onChange={(cancer) => handleCancerChange(cancer)}
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
                            calculateR('cosineSimilarityPublic', {
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
                    <div id="pubPlot">
                      <div style={{ display: pubErr ? 'block' : 'none' }}>
                        <p>
                          An error has occured. Check the debug section for more
                          info.
                        </p>
                      </div>
                      <div style={{ display: pubPlotURL ? 'block' : 'none' }}>
                        <div className="d-flex">
                          <a
                            className="px-2 py-1"
                            href={pubPlotURL}
                            download={pubPlotURL.split('/').slice(-1)[0]}
                          >
                            Download Plot
                          </a>
                          <span className="ml-auto">
                            <Button
                              className="px-2 py-1"
                              variant="link"
                              onClick={() => downloadResults(pubTxtPath)}
                            >
                              Download Results
                            </Button>
                          </span>
                        </div>
                        <div className="p-2 border rounded">
                          <Row>
                            <Col>
                              <img
                                className="w-100 my-4 h-500"
                                src={pubPlotURL}
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
          dispatchCosineSimilarity({
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
