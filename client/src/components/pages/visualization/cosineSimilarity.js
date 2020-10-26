import React from 'react';
import { Form, Row, Col, Button, Accordion, Card } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import {
  dispatchError,
  dispatchCosineSimilarity,
} from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';
import Debug from '../../controls/debug/debug';
import Select from '../../controls/select/select';

const { Group, Label } = Form;
const { Header, Body } = Card;
const { Toggle, Collapse } = Accordion;

export default function CosineSimilarity({ submitR, getRefSigOptions }) {
  const {
    source,
    study,
    studyOptions,
    cancerType,
    pubExperimentalStrategy,
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
  } = useSelector((state) => state.cosineSimilarity);

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
        const response = await fetch(`${rootURL}results/${plotPath}`);
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

  function dispatchRError(fn, status, err = '') {
    if (fn.includes('cosineSimilarityWithin')) {
      dispatchCosineSimilarity({
        withinErr: status,
        debugR: err,
      });
    } else if (fn.includes('cosineSimilarityRefSig')) {
      dispatchCosineSimilarity({
        refErr: status,
        debugR: err,
      });
    } else {
      dispatchCosineSimilarity({
        pubErr: status,
        debugR: err,
      });
    }
  }

  function dispatchOverlay(fn, status) {
    if (fn.includes('cosineSimilarityWithin')) {
      dispatchCosineSimilarity({ withinSubmitOverlay: status });
    } else if (fn.includes('cosineSimilarityRefSig')) {
      dispatchCosineSimilarity({ refSubmitOverlay: status });
    } else {
      dispatchCosineSimilarity({ pubSubmitOverlay: status });
    }
  }

  async function calculateR(fn, args) {
    dispatchOverlay(fn, true);
    dispatchRError(fn, false);

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();
        dispatchOverlay(fn, false);
        dispatchRError(fn, true, err);
      } else {
        const { debugR, output } = await response.json();

        if (Object.keys(output).length) {
          dispatchOverlay(fn, false);
          dispatchRError(fn, false, debugR);
          if (fn.includes('cosineSimilarityWithin')) {
            dispatchCosineSimilarity({
              withinPlotPath: output.plotPath,
              withinTxtPath: output.txtPath,
            });
            setRPlot(output.plotPath, 'within');
          } else if (fn.includes('cosineSimilarityRefSig')) {
            dispatchCosineSimilarity({
              refPlotPath: output.plotPath,
              refTxtPath: output.txtPath,
            });
            setRPlot(output.plotPath, 'refsig');
          } else {
            dispatchCosineSimilarity({
              pubPlotPath: output.plotPath,
              pubTxtPath: output.txtPath,
            });
            setRPlot(output.plotPath, 'pub');
          }
        } else {
          dispatchRError(fn, true, debugR);
          dispatchOverlay(fn, false);
        }
      }
    } catch (err) {
      dispatchError(err);
      dispatchOverlay(fn, false);
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
                      <Select
                        id="csProfileType"
                        label="Profile Type"
                        value={withinProfileType}
                        options={profileOptions}
                        onChange={handleWithinProfileType}
                      />
                    </Col>
                    <Col sm="5">
                      <Select
                        id="csMatrixSize"
                        label="Matrix Size"
                        value={withinMatrixSize}
                        options={withinMatrixOptions}
                        onChange={(matrix) =>
                          dispatchCosineSimilarity({
                            withinMatrixSize: matrix,
                          })
                        }
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
                              experimentalStrategy: pubExperimentalStrategy,
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
                      <Plot
                        plotName={withinPlotPath.split('/').slice(-1)[0]}
                        plotURL={withinPlotURL}
                        txtPath={withinTxtPath}
                      />
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
                      <Select
                        id="csRefProfileType"
                        label="Profile Type"
                        value={refProfileType}
                        options={profileOptions}
                        onChange={(refProfileType) => {
                          dispatchCosineSimilarity({
                            refProfileType: refProfileType,
                          });
                          getSignatureSet(refProfileType);
                        }}
                      />
                    </Col>
                    <Col sm="5">
                      <Select
                        id="csRefSignatureSet"
                        label="Reference Signature Set"
                        value={refSignatureSet}
                        options={refSignatureSetOptions}
                        onChange={(refSignatureSet) => {
                          dispatchCosineSimilarity({
                            refSignatureSet: refSignatureSet,
                          });
                        }}
                      />
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
                              experimentalStrategy: pubExperimentalStrategy,
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
                      <Plot
                        plotName={refPlotPath.split('/').slice(-1)[0]}
                        plotURL={refPlotURL}
                        txtPath={refTxtPath}
                      />
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
                        <Select
                          id="csUserProfileType"
                          label="Profile Type"
                          value={userProfileType}
                          options={profileOptions}
                          onChange={handlePublicProfileType}
                        />
                      </Col>
                      <Col sm="2">
                        <Select
                          id="csUserMatrixSize"
                          label="Matrix Size"
                          value={userMatrixSize}
                          options={userMatrixOptions}
                          onChange={(matrix) =>
                            dispatchCosineSimilarity({
                              userMatrixSize: matrix,
                            })
                          }
                        />
                      </Col>
                      <Col sm="2">
                        <Select
                          id="csPubStudy"
                          label="Study"
                          value={pubStudy}
                          options={studyOptions}
                          onChange={handleStudyChange}
                        />
                      </Col>
                      <Col sm="4">
                        <Select
                          id="csPubCancerType"
                          label="Cancer Type"
                          value={pubCancerType}
                          options={pubCancerTypeOptions}
                          onChange={handleCancerChange}
                        />
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
                        <Plot
                          plotName={pubPlotURL.split('/').slice(-1)[0]}
                          plotURL={pubPlotURL}
                          txtPath={pubTxtPath}
                        />
                      </div>
                    </div>
                  </div>
                </Form>
              </Body>
            </Collapse>
          </Card>
        </Accordion>
      )}
      <Debug msg={debugR} />
    </div>
  );
}
