import React, { useEffect } from 'react';
import {
  Form,
  Row,
  Col,
  Button,
  Popover,
  OverlayTrigger,
  Accordion,
  Card,
} from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import {
  faInfoCircle,
  faPlus,
  faMinus,
} from '@fortawesome/free-solid-svg-icons';
import {
  dispatchError,
  dispatchProfileComparison,
} from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';
import Debug from '../../controls/debug/debug';
import Select from '../../controls/select/select';

const { Group, Label, Control, Text } = Form;
const { Title, Content } = Popover;
const { Header, Body } = Card;
const { Toggle, Collapse } = Accordion;

export default function ProfileComparison({ submitR, getRefSigOptions }) {
  const {
    source,
    study,
    studyOptions,
    cancerType,
    pubExperimentalStrategy,
    pDataOptions,
  } = useSelector((state) => state.visualize);
  const { matrixList, svgList } = useSelector(
    (state) => state.visualizeResults
  );
  const { profileOptions } = useSelector((state) => state.mutationalProfiles);

  const {
    withinProfileType,
    withinSampleName1,
    withinSampleName2,
    sampleOptions,
    refProfileType,
    refSampleName,
    refSignatureSet,
    refSignatureSetOptions,
    refSignatures,
    refCompare,
    userProfileType,
    userSampleName,
    userMatrixSize,
    userMatrixOptions,
    pubSampleName,
    pubSampleOptions,
    pubStudy,
    pubCancerType,
    pubCancerTypeOptions,
    withinPlotPath,
    refPlotPath,
    pubPlotPath,
    withinPlotURL,
    refPlotURL,
    pubPlotURL,
    pubTxtURL,
    displayWithin,
    displayRefSig,
    displayPublic,
    withinErr,
    refErr,
    pubErr,
    debugR,

    withinSubmitOverlay,
    refSubmitOverlay,
    pubSubmitOverlay,
  } = useSelector((state) => state.profileComparison);

  const popover = (
    <Popover id="popover-basic">
      <Title as="h3">{refSignatureSet}</Title>
      <Content>
        {refSignatures.length > 0 ? (
          <ul>
            {refSignatures.map((signature) => (
              <li key={signature}>{signature}</li>
            ))}
          </ul>
        ) : (
          <p>No Signature Available</p>
        )}
      </Content>
    </Popover>
  );

  // get inital signatures from initially selected signature set
  useEffect(() => {
    if (refProfileType && refSignatureSet && !refSignatures.length) {
      getSignatures(refProfileType, refSignatureSet);
    }
  }, [refProfileType, refSignatureSet]);

  // get public data samples
  useEffect(() => {
    if (!pubSampleOptions.length && source == 'user' && pubStudy) {
      getPublicSamples(pubStudy, pubCancerType);
    }
  }, [pDataOptions, pubStudy]);

  function handleOverlay(fn, status) {
    if (fn.includes('profileComparisonWithin')) {
      dispatchProfileComparison({ withinSubmitOverlay: status });
    } else if (fn.includes('profileComparisonRefSig')) {
      dispatchProfileComparison({ refSubmitOverlay: status });
    } else {
      dispatchProfileComparison({ pubSubmitOverlay: status });
    }
  }

  async function setRPlot(plotPath, fn) {
    handleOverlay(fn, true);
    if (plotPath) {
      try {
        const response = await fetch(`api/results/${plotPath}`);
        if (!response.ok) {
          // console.log(await response.json());
        } else {
          const pic = await response.blob();
          const objectURL = URL.createObjectURL(pic);

          if (fn.includes('profileComparisonWithin')) {
            if (withinPlotURL) URL.revokeObjectURL(withinPlotURL);
            dispatchProfileComparison({
              withinPlotURL: objectURL,
            });
          } else if (fn.includes('profileComparisonRefSig')) {
            if (refPlotURL) URL.revokeObjectURL(refPlotURL);
            dispatchProfileComparison({
              refPlotURL: objectURL,
            });
          } else {
            if (pubPlotURL) URL.revokeObjectURL(pubPlotURL);
            dispatchProfileComparison({
              pubPlotURL: objectURL,
            });
          }
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (fn.includes('profileComparisonWithin')) {
        if (withinPlotURL) URL.revokeObjectURL(withinPlotURL);
        dispatchProfileComparison({ withinErr: true, withinPlotURL: '' });
      } else if (fn.includes('profileComparisonRefSig')) {
        if (refPlotURL) URL.revokeObjectURL(refPlotURL);
        dispatchProfileComparison({ refErr: true, refPlotURL: '' });
      } else {
        if (pubPlotURL) URL.revokeObjectURL(pubPlotURL);
        dispatchProfileComparison({ pubErr: true, pubPlotURL: '' });
      }
    }
    handleOverlay(fn, false);
  }

  // get Signature Reference Sets for dropdown options
  async function getSignatureSet(profileType) {
    if (profileType) {
      dispatchProfileComparison({ refSubmitOverlay: true });
      try {
        const response = await getRefSigOptions(profileType);

        if (response.ok) {
          const signatureSetOptions = await response.json();

          dispatchProfileComparison({
            refSignatureSetOptions: signatureSetOptions,
            refSignatureSet: signatureSetOptions[0],
            refSubmitOverlay: false,
          });
          getSignatures(profileType, signatureSetOptions[0]);
        } else {
          dispatchError(await response.json());
          dispatchProfileComparison({ refSubmitOverlay: false });
        }
      } catch (err) {
        dispatchError(err);
        dispatchProfileComparison({ refSubmitOverlay: false });
      }
    }
  }

  // get signature options for compare
  async function getSignatures(profileType, signatureSetName) {
    if (signatureSetName) {
      dispatchProfileComparison({ refSubmitOverlay: true });
      try {
        const response = await fetch(`api/visualizeR/getSignatures`, {
          method: 'POST',
          headers: {
            Accept: 'application/json',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            profileType: profileType,
            signatureSetName: signatureSetName,
          }),
        });

        if (response.ok) {
          const signatures = await response.json();

          dispatchProfileComparison({
            refSignatures: signatures,
            refCompare: signatures[0],
            refSubmitOverlay: false,
          });
        } else {
          dispatchError(await response.json());
          dispatchProfileComparison({ refSubmitOverlay: false });
        }
      } catch (err) {
        dispatchError(err);
        dispatchProfileComparison({ refSubmitOverlay: false });
      }
    }
  }

  async function calculateR(fn, args) {
    handleOverlay(fn, true);
    if (fn.includes('profileComparisonWithin')) {
      dispatchProfileComparison({
        withinErr: false,
        debugR: '',
      });
    } else if (fn.includes('profileComparisonRefSig')) {
      dispatchProfileComparison({
        refErr: false,
        debugR: '',
      });
    } else {
      dispatchProfileComparison({
        pubErr: false,
        debugR: '',
      });
    }
    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();
        dispatchProfileComparison({ debugR: err });
        handleOverlay(fn, false);
      } else {
        const { debugR, output } = await response.json();

        dispatchProfileComparison({ debugR: debugR });
        if (Object.keys(output).length) {
          if (fn.includes('profileComparisonWithin')) {
            dispatchProfileComparison({ withinPlotPath: output.plotPath });
            setRPlot(output.plotPath, fn);
          } else if (fn.includes('profileComparisonRefSig')) {
            dispatchProfileComparison({ refPlotPath: output.plotPath });
            setRPlot(output.plotPath, fn);
          } else {
            dispatchProfileComparison({ pubPlotPath: output.plotPath });
            setRPlot(output.plotPath, fn);
          }
        } else {
          handleOverlay(fn, false);
          if (fn.includes('profileComparisonWithin')) {
            dispatchProfileComparison({ withinPlotPath: '', withinErr: true });
          } else if (fn.includes('profileComparisonRefSig')) {
            dispatchProfileComparison({ refPlotPath: '', refErr: true });
          } else {
            dispatchProfileComparison({ pubPlotPath: '', pubErr: true });
          }
        }
      }
    } catch (err) {
      dispatchError(err);
      handleOverlay(fn, false);
    }
  }

  async function getPublicSamples(study, cancerType) {
    const args = {
      study: study,
      cancerType: cancerType,
      experimentalStrategy: [
        ...new Set(
          pDataOptions
            .filter(
              (data) => data.Study == study && data.Cancer_Type == cancerType
            )
            .map((data) => data.Dataset)
        ),
      ][0],
    };

    handleOverlay('pub', true);
    try {
      const response = await fetch(`api/getPublicData`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(args),
      });

      if (response.ok) {
        const { svgList } = await response.json();
        const pubSamples = [...new Set(svgList.map((plot) => plot.Sample))];

        dispatchProfileComparison({
          pubSampleOptions: pubSamples,
          pubSampleName: pubSamples[0],
        });
      } else {
        dispatchProfileComparison({
          pubSampleOptions: ['NA'],
          pubErr: 'Error retrieving pub data samples',
        });
      }
    } catch (err) {
      dispatchError(err);
      dispatchProfileComparison({
        pubErr: 'Error retrieving pub data samples',
      });
    }
    handleOverlay('pub', false);
  }

  function handleStudyChange(study) {
    const cancerTypeOptions = [
      ...new Set(
        pDataOptions
          .filter((data) => data.Study == study)
          .map((data) => data.Cancer_Type)
      ),
    ];

    dispatchProfileComparison({
      pubStudy: study,
      pubCancerType: cancerTypeOptions[0],
      pubCancerTypeOptions: cancerTypeOptions,
    });
    getPublicSamples(study, cancerTypeOptions[0]);
  }

  function handleCancerChange(cancer) {
    dispatchProfileComparison({
      pubCancerType: cancer,
    });
    getPublicSamples(pubStudy, cancer);
  }

  function handleProfile(profile) {
    const matrixOptions = [
      ...new Set(
        svgList
          .filter((plot) => plot.Profile_Type == profile)
          .map((plot) => plot.Matrix_Size)
      ),
    ];

    dispatchProfileComparison({
      userProfileType: profile,
      userMatrixOptions: matrixOptions,
      userMatrixSize: matrixOptions[0],
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
              dispatchProfileComparison({
                displayWithin: !displayWithin,
              })
            }
          >
            {displayWithin == true ? (
              <FontAwesomeIcon icon={faMinus} />
            ) : (
              <FontAwesomeIcon icon={faPlus} />
            )}{' '}
            Comparison Within Samples
          </Toggle>
          <Collapse eventKey="0">
            <Body>
              <Form>
                <LoadingOverlay active={withinSubmitOverlay} />
                <Row className="justify-content-center">
                  <Col sm="1">
                    <Select
                      disabled={sampleOptions.length < 2}
                      id="pcProfileTypeWithin"
                      label="Profile Type"
                      value={withinProfileType}
                      options={profileOptions}
                      onChange={(profile) =>
                        dispatchProfileComparison({
                          withinProfileType: profile,
                        })
                      }
                    />
                  </Col>
                  <Col sm="5">
                    <Select
                      disabled={sampleOptions.length < 2}
                      id="pcSample1"
                      label="Sample Name 1"
                      value={withinSampleName1}
                      options={sampleOptions}
                      onChange={(name) => {
                        dispatchProfileComparison({
                          withinSampleName1: name,
                        });
                      }}
                    />
                  </Col>
                  <Col sm="5">
                    <Select
                      disabled={sampleOptions.length < 2}
                      id="pcSample2"
                      label="Sample Name 2"
                      value={withinSampleName2}
                      options={sampleOptions}
                      onChange={(name) => {
                        dispatchProfileComparison({
                          withinSampleName2: name,
                        });
                      }}
                    />
                  </Col>
                  <Col sm="1" className="m-auto">
                    <Button
                      disabled={sampleOptions.length < 2}
                      variant="primary"
                      onClick={() => {
                        if (source == 'user') {
                          calculateR('profileComparisonWithin', {
                            profileType: withinProfileType,
                            sampleName1: withinSampleName1,
                            sampleName2: withinSampleName2,
                            matrixList: JSON.stringify(
                              matrixList.filter(
                                (matrix) =>
                                  matrix.Profile_Type == withinProfileType
                              )
                            ),
                          });
                        } else {
                          calculateR('profileComparisonWithinPublic', {
                            profileType: withinProfileType,
                            sampleName1: withinSampleName1,
                            sampleName2: withinSampleName2,
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
                {sampleOptions.length < 2 && (
                  <Row>
                    <Col>Unavailable - More than one Sample Required</Col>
                  </Row>
                )}

                <div id="pcWithinPlot">
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
                    />
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
              dispatchProfileComparison({
                displayRefSig: !displayRefSig,
              })
            }
          >
            {displayRefSig == true ? (
              <FontAwesomeIcon icon={faMinus} />
            ) : (
              <FontAwesomeIcon icon={faPlus} />
            )}{' '}
            Comparison to Reference Signatures
          </Toggle>
          <Collapse eventKey="1">
            <Body>
              <Form className="my-2">
                <LoadingOverlay active={refSubmitOverlay} />
                <div>
                  <Row className="justify-content-center">
                    <Col sm="1">
                      <Select
                        id="pcProfileTypeRef"
                        label="Profile Type"
                        value={refProfileType}
                        options={profileOptions}
                        onChange={(refProfileType) => {
                          dispatchProfileComparison({
                            refProfileType: refProfileType,
                          });
                          getSignatureSet(refProfileType);
                        }}
                      />
                    </Col>
                    <Col sm="3">
                      <Select
                        id="sampleNameRefSig"
                        label="Sample Name"
                        value={refSampleName}
                        options={sampleOptions}
                        onChange={(name) => {
                          dispatchProfileComparison({
                            refSampleName: name,
                          });
                        }}
                      />
                    </Col>
                    <Col sm="4">
                      <Select
                        id="pcRefSet"
                        label="Reference Signature Set"
                        value={refSignatureSet}
                        options={refSignatureSetOptions}
                        onChange={(refSignatureSet) => {
                          dispatchProfileComparison({
                            refSignatureSet: refSignatureSet,
                          });
                          getSignatures(refProfileType, refSignatureSet);
                        }}
                      />
                    </Col>
                    <Col sm="3">
                      <Group controlId="signatureSet">
                        <Label>
                          Compare Signatures{' '}
                          <OverlayTrigger
                            trigger="click"
                            placement="top"
                            overlay={popover}
                            rootClose
                          >
                            <Button
                              aria-label="compare signatures info"
                              variant="link"
                              className="p-0 font-weight-bold "
                            >
                              <FontAwesomeIcon
                                icon={faInfoCircle}
                                style={{ verticalAlign: 'baseline' }}
                              />
                            </Button>
                          </OverlayTrigger>
                        </Label>
                        <Control
                          value={refCompare}
                          onChange={(e) => {
                            dispatchProfileComparison({
                              refCompare: e.target.value,
                            });
                          }}
                        ></Control>
                        <Text className="text-muted">
                          (Ex. 0.8*SBS5;0.1*SBS1)
                        </Text>
                      </Group>
                    </Col>
                    <Col sm="1" className="m-auto">
                      <Button
                        variant="primary"
                        onClick={() => {
                          if (source == 'user') {
                            calculateR('profileComparisonRefSig', {
                              profileType: refProfileType,
                              sampleName: refSampleName,
                              signatureSet: refSignatureSet,
                              compare: refCompare,
                              matrixList: JSON.stringify(
                                matrixList.filter(
                                  (matrix) =>
                                    matrix.Profile_Type == withinProfileType
                                )
                              ),
                            });
                          } else {
                            calculateR('profileComparisonRefSigPublic', {
                              profileType: refProfileType,
                              sampleName: refSampleName,
                              signatureSet: refSignatureSet,
                              compare: refCompare,
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

                  <div id="refPlotDownload">
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
                dispatchProfileComparison({
                  displayPublic: !displayPublic,
                })
              }
            >
              {displayPublic == true ? (
                <FontAwesomeIcon icon={faMinus} />
              ) : (
                <FontAwesomeIcon icon={faPlus} />
              )}{' '}
              Comparison to Public Data
            </Toggle>
            <Collapse eventKey="2">
              <Body>
                <Form>
                  <LoadingOverlay active={pubSubmitOverlay} />
                  <div>
                    <Row className="justify-content-center">
                      <Col sm="1">
                        <Select
                          id="pcUserProfileType"
                          label="Profile Type"
                          value={userProfileType}
                          options={profileOptions}
                          onChange={handleProfile}
                        />
                      </Col>
                      <Col sm="1">
                        <Select
                          id="pcUserMatrixSize"
                          label="Matrix Size"
                          value={userMatrixSize}
                          options={userMatrixOptions}
                          onChange={(matrix) =>
                            dispatchProfileComparison({
                              userMatrixSize: matrix,
                            })
                          }
                        />
                      </Col>
                      <Col sm="2">
                        <Select
                          id="pcUserSampleName"
                          label="Sample Name"
                          value={userSampleName}
                          options={sampleOptions}
                          onChange={(name) => {
                            dispatchProfileComparison({
                              userSampleName: name,
                            });
                          }}
                        />
                      </Col>
                      <Col sm="2">
                        <Select
                          id="pcPubStudy"
                          label="Study"
                          value={pubStudy}
                          options={studyOptions}
                          onChange={handleStudyChange}
                        />
                      </Col>
                      <Col sm="2">
                        <Select
                          id="pcPubCancerType"
                          label="Cancer Type"
                          value={pubCancerType}
                          options={pubCancerTypeOptions}
                          onChange={handleCancerChange}
                        />
                      </Col>
                      <Col sm="3">
                        <Select
                          id="pcPubSampleName"
                          label="Public Sample Name"
                          value={pubSampleName}
                          options={pubSampleOptions}
                          onChange={(name) => {
                            dispatchProfileComparison({
                              pubSampleName: name,
                            });
                          }}
                        />
                      </Col>
                      <Col sm="1" className="m-auto">
                        <Button
                          variant="primary"
                          onClick={() =>
                            calculateR('profileComparisonPublic', {
                              profileName: userProfileType + userMatrixSize,
                              matrixFile: matrixList.filter(
                                (path) =>
                                  path.Profile_Type == userProfileType &&
                                  path.Matrix_Size == userMatrixSize
                              )[0].Path,
                              userSample: userSampleName,
                              study: pubStudy,
                              cancerType: pubCancerType,
                              publicSample: pubSampleName,
                            })
                          }
                        >
                          Calculate
                        </Button>
                      </Col>
                    </Row>

                    <div id="pcPubPlot">
                      <div style={{ display: pubErr ? 'block' : 'none' }}>
                        <p>
                          An error has occured. Check the debug section for more
                          info.
                        </p>
                      </div>
                      <div style={{ display: pubPlotURL ? 'block' : 'none' }}>
                        <Plot
                          plotName={pubPlotPath.split('/').slice(-1)[0]}
                          plotURL={pubPlotURL}
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
