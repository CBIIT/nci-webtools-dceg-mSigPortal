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
  dispatchMutationalPattern,
} from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

const { Group, Label, Control, Text } = Form;
const { Header, Body } = Card;
const { Toggle, Collapse } = Accordion;

export default function MutationalPattern({
  downloadResults,
  submitR,
  getRefSigOptions,
}) {
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
    profileType,
    matrixSize,
    matrixOptions,
    proportion,
    pattern,
    pubStudy,
    pubCancerType,
    pubCancerTypeOptions,
    plotPath,

    plotURL,
    display,
    err,
    debugR,
    displayDebug,
    submitOverlay,
  } = useSelector((state) => state.mutationalPattern);

  const selectFix = {
    styles: {
      menuPortal: (base) => ({ ...base, zIndex: 9999 }),
    },
    menuPortalTarget: document.body,
  };

  function setOverlay(type, display) {
    if (type == 'within') {
      dispatchMutationalPattern({ submitOverlay: display });
    } else if (type == 'refsig') {
      dispatchMutationalPattern({ refSubmitOverlay: display });
    } else {
      dispatchMutationalPattern({ pubSubmitOverlay: display });
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

          if (plotURL) URL.revokeObjectURL(plotURL);
          dispatchMutationalPattern({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchMutationalPattern({ err: true, plotURL: '' });
    }
    setOverlay(type, false);
  }

  // get Signature Reference Sets for dropdown options
  async function getSignatureSet(profileType) {
    if (profileType) {
      dispatchMutationalPattern({ refSubmitOverlay: true });
      try {
        const response = await getRefSigOptions(profileType);

        if (response.ok) {
          const signatureSetOptions = await response.json();

          dispatchMutationalPattern({
            refSignatureSetOptions: signatureSetOptions,
            refSignatureSet: signatureSetOptions[0],
            refSubmitOverlay: false,
          });
          getSignatures(profileType, signatureSetOptions[0]);
        } else {
          dispatchError(await response.json());
          dispatchMutationalPattern({ refSubmitOverlay: false });
        }
      } catch (err) {
        dispatchError(err);
        dispatchMutationalPattern({ refSubmitOverlay: false });
      }
    }
  }

  // get signature options for compare
  async function getSignatures(profileType, signatureSetName) {
    if (signatureSetName) {
      dispatchMutationalPattern({ refSubmitOverlay: true });
      try {
        const response = await fetch(`${rootURL}api/visualizeR/getSignatures`, {
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

          dispatchMutationalPattern({
            refSignatures: signatures,
            refCompare: signatures[0],
            refSubmitOverlay: false,
          });
        } else {
          dispatchError(await response.json());
          dispatchMutationalPattern({ refSubmitOverlay: false });
        }
      } catch (err) {
        dispatchError(err);
        dispatchMutationalPattern({ refSubmitOverlay: false });
      }
    }
  }

  async function calculateR(fn, args) {
    dispatchMutationalPattern({
      submitOverlay: true,
      err: false,
      debugR: '',
    });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();
        dispatchMutationalPattern({ debugR: err });

        dispatchMutationalPattern({ submitOverlay: false });
      } else {
        const { debugR, output } = await response.json();

        dispatchMutationalPattern({ debugR: debugR });

        dispatchMutationalPattern({ plotPath: output.plotPath });
        setRPlot(output.plotPath, 'within');
      }
    } catch (err) {
      dispatchError(err);

      dispatchMutationalPattern({ submitOverlay: false });
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
      const response = await fetch(`${rootURL}visualize/getPublicData`, {
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

        dispatchMutationalPattern({
          pubSampleOptions: pubSamples,
          pubSampleName: pubSamples[0],
        });
      } else {
        dispatchMutationalPattern({
          pubSampleOptions: ['NA'],
          pubErr: 'Error retrieving pub data samples',
        });
      }
    } catch (err) {
      dispatchError(err);
      dispatchMutationalPattern({
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

    dispatchMutationalPattern({
      pubStudy: study,
      pubCancerType: cancerTypeOptions[0],
      pubCancerTypeOptions: cancerTypeOptions,
    });
    getPublicSamples(study, cancerTypeOptions[0]);
  }

  function handleCancerChange(cancer) {
    dispatchMutationalPattern({
      pubCancerType: cancer,
    });
    getPublicSamples(pubStudy, cancer);
  }

  function handleProfileType(profileType) {
    const matrixOptions = [
      ...new Set(
        matrixList
          .filter((matrix) => matrix.Profile_Type == profileType)
          .map((matrix) => matrix.Matrix_Size)
      ),
    ];

    dispatchMutationalPattern({
      profileType: profileType,
      matrixSize: matrixOptions[0],
      matrixOptions: matrixOptions,
    });
  }

  return (
    <div>
      {source == 'user' && (
        <Accordion defaultActiveKey="0">
          <Card>
            <Toggle
              className="font-weight-bold"
              as={Header}
              eventKey="0"
              onClick={() =>
                dispatchMutationalPattern({
                  display: !display,
                })
              }
            >
              {display == true ? (
                <FontAwesomeIcon icon={faMinus} />
              ) : (
                <FontAwesomeIcon icon={faPlus} />
              )}{' '}
              Motif Analysis
            </Toggle>
            <Collapse eventKey="0">
              <Body>
                <Form>
                  <LoadingOverlay active={submitOverlay} />
                  <div>
                    <Row className="justify-content-center">
                      <Col sm="2">
                        <Group controlId="profileType">
                          <Label>Profile Type</Label>
                          <Select
                            options={profileOptions}
                            value={[profileType]}
                            onChange={(profile) => handleProfileType(profile)}
                            getOptionLabel={(option) => option}
                            getOptionValue={(option) => option}
                            {...selectFix}
                          />
                        </Group>
                      </Col>
                      <Col sm="2">
                        <Label>Matrix Size</Label>
                        <Select
                          options={matrixOptions}
                          value={[matrixSize]}
                          onChange={(matrix) =>
                            dispatchMutationalPattern({
                              matrixSize: matrix,
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
                      <Col sm="2">
                        <Label>
                          Minimal Proportion mutations within Each Mutational
                          Pattern
                        </Label>
                        <Control
                          value={proportion}
                          onChange={(e) => {
                            dispatchMutationalPattern({
                              proportion: e.target.value,
                            });
                          }}
                        ></Control>{' '}
                        <Text className="text-muted">(Ex. 0.8)</Text>
                      </Col>
                      <Col sm="2">
                        <Label>Mutational Pattern</Label>
                        <Control
                          value={pattern}
                          onChange={(e) => {
                            dispatchMutationalPattern({
                              pattern: e.target.value,
                            });
                          }}
                        ></Control>
                        <Text className="text-muted">(Ex. NCG>NTG)</Text>
                      </Col>

                      <Col sm="2" className="m-auto">
                        <Button
                          variant="primary"
                          onClick={() => {
                            if (source == 'user') {
                              calculateR('mutationalPattern', {
                                matrixFile: matrixList.filter(
                                  (path) =>
                                    path.Profile_Type == profileType &&
                                    path.Matrix_Size == matrixSize
                                )[0].Path,
                                study: pubStudy,
                                proportion: proportion,
                                pattern: pattern,
                              });
                            } else {
                              calculateR('profileComparisonWithinPublic', {
                                profileType: profileType,
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

                    <div id="plot">
                      <div style={{ display: err ? 'block' : 'none' }}>
                        <p>
                          An error has occured. Check the debug section for more
                          info.
                        </p>
                      </div>
                      <div style={{ display: plotURL ? 'block' : 'none' }}>
                        <div className="d-flex">
                          <a
                            className="px-2 py-1"
                            href={plotURL}
                            download={plotURL.split('/').slice(-1)[0]}
                          >
                            Download Plot
                          </a>
                        </div>
                        <div className="p-2 border rounded">
                          <Row>
                            <Col>
                              <img className="w-100 my-4" src={plotURL}></img>
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

      {source == 'public' && (
        <Accordion defaultActiveKey="1">
          <Card>
            <Toggle
              className="font-weight-bold"
              as={Header}
              eventKey="1"
              onClick={() =>
                dispatchMutationalPattern({
                  display: !display,
                })
              }
            >
              {display == true ? (
                <FontAwesomeIcon icon={faMinus} />
              ) : (
                <FontAwesomeIcon icon={faPlus} />
              )}{' '}
              Motif Analysis
            </Toggle>
            <Collapse eventKey="1">
              <Body>
                <Form>
                  <LoadingOverlay active={submitOverlay} />
                  <div>
                    <Row className="justify-content-center">
                      <Col sm="2">
                        <Group controlId="profileType">
                          <Label>Profile Type</Label>
                          <Select
                            options={profileOptions}
                            value={[profileType]}
                            onChange={(profile) => handleProfileType(profile)}
                            getOptionLabel={(option) => option}
                            getOptionValue={(option) => option}
                            {...selectFix}
                          />
                        </Group>
                      </Col>
                      <Col sm="2">
                        <Label>Matrix Size</Label>
                        <Select
                          options={matrixOptions}
                          value={[matrixSize]}
                          onChange={(matrix) =>
                            dispatchMutationalPattern({
                              matrixSize: matrix,
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
                      <Col sm="2">
                        <Label>
                          Minimal Proportion mutations within Each Mutational
                          Pattern
                        </Label>
                        <Control
                          value={proportion}
                          onChange={(e) => {
                            dispatchMutationalPattern({
                              proportion: e.target.value,
                            });
                          }}
                        ></Control>{' '}
                        <Text className="text-muted">(Ex. 0.8)</Text>
                      </Col>
                      <Col sm="2">
                        <Label>Mutational Pattern</Label>
                        <Control
                          value={pattern}
                          onChange={(e) => {
                            dispatchMutationalPattern({
                              pattern: e.target.value,
                            });
                          }}
                        ></Control>
                        <Text className="text-muted">(Ex. NCG>NTG)</Text>
                      </Col>

                      <Col sm="2" className="m-auto">
                        <Button
                          variant="primary"
                          onClick={() => {
                            if (source == 'user') {
                              calculateR('mutationalPattern', {
                                matrixFile: matrixList.filter(
                                  (path) =>
                                    path.Profile_Type == profileType &&
                                    path.Matrix_Size == matrixSize
                                )[0].Path,
                                study: pubStudy,
                                proportion: proportion,
                                pattern: pattern,
                              });
                            } else {
                              calculateR('profileComparisonWithinPublic', {
                                profileType: profileType,
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

                    <div id="plot">
                      <div style={{ display: err ? 'block' : 'none' }}>
                        <p>
                          An error has occured. Check the debug section for more
                          info.
                        </p>
                      </div>
                      <div style={{ display: plotURL ? 'block' : 'none' }}>
                        <div className="d-flex">
                          <a
                            className="px-2 py-1"
                            href={plotURL}
                            download={plotURL.split('/').slice(-1)[0]}
                          >
                            Download Plot
                          </a>
                        </div>
                        <div className="p-2 border rounded">
                          <Row>
                            <Col>
                              <img className="w-100 my-4" src={plotURL}></img>
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
          dispatchMutationalPattern({
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
