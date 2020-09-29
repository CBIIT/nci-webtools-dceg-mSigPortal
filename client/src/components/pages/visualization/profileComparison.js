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
import Select from 'react-select';
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
  const { matrixList } = useSelector((state) => state.visualizeResults);
  const { nameOptions, profileOptions } = useSelector(
    (state) => state.mutationalProfiles
  );
  const rootURL = window.location.pathname;
  const {
    withinProfileType,
    withinSampleName1,
    withinSampleName2,
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
    displayDebug,
    withinSubmitOverlay,
    refSubmitOverlay,
    pubSubmitOverlay,
  } = useSelector((state) => state.profileComparison);

  const selectFix = {
    styles: {
      menuPortal: (base) => ({ ...base, zIndex: 9999 }),
    },
    menuPortalTarget: document.body,
  };

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

  function setOverlay(type, display) {
    if (type == 'within') {
      dispatchProfileComparison({ withinSubmitOverlay: display });
    } else if (type == 'refsig') {
      dispatchProfileComparison({ refSubmitOverlay: display });
    } else {
      dispatchProfileComparison({ pubSubmitOverlay: display });
    }
  }

  async function setRPlot(plotPath, type) {
    setOverlay(type, true);
    if (plotPath) {
      try {
        const response = await fetch(`${rootURL}getSVG`, {
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
            dispatchProfileComparison({
              withinPlotURL: objectURL,
            });
          } else if (type == 'refsig') {
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
      if (type == 'within') {
        if (withinPlotURL) URL.revokeObjectURL(withinPlotURL);
        dispatchProfileComparison({ withinErr: true, withinPlotURL: '' });
      } else if (type == 'refsig') {
        if (refPlotURL) URL.revokeObjectURL(refPlotURL);
        dispatchProfileComparison({ refErr: true, refPlotURL: '' });
      } else {
        if (pubPlotURL) URL.revokeObjectURL(pubPlotURL);
        dispatchProfileComparison({ pubErr: true, pubPlotURL: '' });
      }
    }
    setOverlay(type, false);
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
        const response = await fetch(`${rootURL}visualizeR/getSignatures`, {
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
    if (fn.includes('profileComparisonWithin')) {
      dispatchProfileComparison({
        withinSubmitOverlay: true,
        withinErr: false,
        debugR: '',
      });
    } else if (fn.includes('profileComparisonRefSig')) {
      dispatchProfileComparison({
        refSubmitOverlay: true,
        refErr: false,
        debugR: '',
      });
    } else {
      dispatchProfileComparison({
        pubSubmitOverlay: true,
        pubErr: false,
        debugR: '',
      });
    }
    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();
        dispatchProfileComparison({ debugR: err });
        if (fn.includes('profileComparisonWithin')) {
          dispatchProfileComparison({ withinSubmitOverlay: false });
        } else if (fn.includes('profileComparisonRefSig')) {
          dispatchProfileComparison({ refSubmitOverlay: false });
        } else {
          dispatchProfileComparison({ pubSubmitOverlay: false });
        }
      } else {
        const { debugR, output } = await response.json();

        dispatchProfileComparison({ debugR: debugR });

        if (fn.includes('profileComparisonWithin')) {
          dispatchProfileComparison({ withinPlotPath: output.plotPath });
          setRPlot(output.plotPath, 'within');
        } else if (fn.includes('profileComparisonRefSig')) {
          dispatchProfileComparison({ refPlotPath: output.plotPath });
          setRPlot(output.plotPath, 'refsig');
        } else {
          dispatchProfileComparison({ pubPlotPath: output.plotPath });
          setRPlot(output.plotPath, 'pub');
        }
      }
    } catch (err) {
      dispatchError(err);
      if (fn.includes('profileComparisonWithin')) {
        dispatchProfileComparison({ withinSubmitOverlay: false });
      } else if (fn.includes('profileComparisonRefSig')) {
        dispatchProfileComparison({ refSubmitOverlay: false });
      } else {
        dispatchProfileComparison({ pubSubmitOverlay: false });
      }
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

    setOverlay('pub', true);
    try {
      const response = await fetch(`${rootURL}getPublicData`, {
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
    setOverlay('pub', false);
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
                <div>
                  <Row className="justify-content-center">
                    <Col sm="1">
                      <Group controlId="profileTypeWithin">
                        <Label>Profile Type</Label>
                        <Select
                          options={profileOptions}
                          value={[withinProfileType]}
                          onChange={(profile) =>
                            dispatchProfileComparison({
                              withinProfileType: profile,
                            })
                          }
                          getOptionLabel={(option) => option}
                          getOptionValue={(option) => option}
                          {...selectFix}
                        />
                      </Group>
                    </Col>
                    <Col sm="5">
                      <Label>Sample Name 1</Label>
                      <Select
                        options={nameOptions}
                        value={[withinSampleName1]}
                        onChange={(name) => {
                          dispatchProfileComparison({
                            withinSampleName1: name,
                          });
                        }}
                        getOptionLabel={(option) => option}
                        getOptionValue={(option) => option}
                        {...selectFix}
                      />
                    </Col>
                    <Col sm="5">
                      <Label>Sample Name 2</Label>
                      <Select
                        options={nameOptions}
                        value={[withinSampleName2]}
                        onChange={(name) => {
                          dispatchProfileComparison({
                            withinSampleName2: name,
                          });
                        }}
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
                      <Group controlId="profileTypeRefSig">
                        <Label>Profile Type</Label>
                        <Select
                          options={profileOptions}
                          value={[refProfileType]}
                          onChange={(refProfileType) => {
                            dispatchProfileComparison({
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
                    <Col sm="3">
                      <Group controlId="sampleNameRefSig">
                        <Label>Sample Name</Label>
                        <Select
                          options={nameOptions}
                          value={[refSampleName]}
                          onChange={(name) => {
                            dispatchProfileComparison({
                              refSampleName: name,
                            });
                          }}
                          getOptionLabel={(option) => option}
                          getOptionValue={(option) => option}
                          {...selectFix}
                        />
                      </Group>
                    </Col>
                    <Col sm="4">
                      <Group controlId="signatureSet">
                        <Label>Reference Signature Set</Label>
                        <Select
                          options={refSignatureSetOptions}
                          value={[refSignatureSet]}
                          onChange={(refSignatureSet) => {
                            dispatchProfileComparison({
                              refSignatureSet: refSignatureSet,
                            });
                            getSignatures(refProfileType, refSignatureSet);
                          }}
                          getOptionLabel={(option) => option}
                          getOptionValue={(option) => option}
                          {...selectFix}
                        />
                      </Group>
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
                        <Group controlId="userProfileType">
                          <Label>Profile Type</Label>
                          <Select
                            options={profileOptions}
                            value={[userProfileType]}
                            onChange={(profile) =>
                              dispatchProfileComparison({
                                userProfileType: profile,
                              })
                            }
                            getOptionLabel={(option) => option}
                            getOptionValue={(option) => option}
                            {...selectFix}
                          />
                        </Group>
                      </Col>
                      <Col sm="1">
                        <Label>Matrix Size</Label>
                        <Select
                          options={userMatrixOptions}
                          value={[userMatrixSize]}
                          onChange={(matrix) =>
                            dispatchProfileComparison({
                              userMatrixSize: matrix,
                            })
                          }
                          getOptionLabel={(option) => option}
                          getOptionValue={(option) => option}
                          {...selectFix}
                        />
                      </Col>
                      <Col sm="2">
                        <Label>Sample Name</Label>
                        <Select
                          options={nameOptions}
                          value={[userSampleName]}
                          onChange={(name) => {
                            dispatchProfileComparison({
                              userSampleName: name,
                            });
                          }}
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
                      <Col sm="2">
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
                      <Col sm="3">
                        <Label>Public Sample Name</Label>
                        <Select
                          options={pubSampleOptions}
                          value={[pubSampleName]}
                          onChange={(name) => {
                            dispatchProfileComparison({
                              pubSampleName: name,
                            });
                          }}
                          getOptionLabel={(option) => option}
                          getOptionValue={(option) => option}
                          {...selectFix}
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

                    <div id="pcWithinPlot">
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

      <Button
        variant="link"
        className="p-0 mt-5"
        onClick={() =>
          dispatchProfileComparison({
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
