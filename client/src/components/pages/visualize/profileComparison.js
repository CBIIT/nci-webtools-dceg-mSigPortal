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

const { Group, Label, Control, Text } = Form;
const { Title, Content } = Popover;
const { Header, Body } = Card;
const { Toggle, Collapse } = Accordion;

export default function ProfileComparison({ submitR, getRefSigOptions }) {
  const { displayTab } = useSelector((state) => state.visualizeResults);
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
    withinPlotPath,
    refPlotPath,
    withinPlotURL,
    refPlotURL,
    displayWithin,
    displayRefSig,
    withinErr,
    refErr,
    debugR,
    displayDebug,
    withinSubmitOverlay,
    refSubmitOverlay,
  } = useSelector((state) => state.profileComparison);

  const selectFix = {
    styles: {
      menuPortal: (base) => ({ ...base, zIndex: 9999 }),
    },
    menuPortalTarget: document.body,
  };

  // get inital signatures from initially selected signature set
  useEffect(() => {
    if (refProfileType && refSignatureSet && !refSignatures.length) {
      getSignatures(refProfileType, refSignatureSet);
    }
  }, [refProfileType, refSignatureSet]);

  // calculate r on load
  // useEffect(() => {
  //   if (
  //     withinProfileType.length &&
  //     withinSampleName1.length &&
  //     withinSampleName2.length &&
  //     !withinPlotPath &&
  //     !withinSubmitOverlay &&
  //     displayTab == 'profileComparison'
  //   ) {
  //     calculateR('profileComparisonWithin', {
  //       profileType: withinProfileType,
  //       sampleName1: withinSampleName1,
  //       sampleName2: withinSampleName2,
  //     });
  //   }
  // }, [displayTab]);

  // useEffect(() => {
  //   if (
  //     refProfileType.length &&
  //     refSampleName.length &&
  //     refSignatureSet.length &&
  //     refCompare.length &&
  //     !refPlotPath &&
  //     !refSubmitOverlay &&
  //     displayTab == 'profileComparison'
  //   ) {
  //     calculateR('profileComparisonRefSig', {
  //       profileType: refProfileType,
  //       sampleName: refSampleName,
  //       signatureSet: refSignatureSet,
  //       compare: refCompare,
  //     });
  //   }
  // }, [displayTab]);

  function setOverlay(type, display) {
    if (type == 'within') {
      dispatchProfileComparison({ withinSubmitOverlay: display });
    } else {
      dispatchProfileComparison({ refSubmitOverlay: display });
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
            dispatchProfileComparison({
              withinPlotURL: objectURL,
            });
          } else {
            if (refPlotURL) URL.revokeObjectURL(refPlotURL);
            dispatchProfileComparison({
              refPlotURL: objectURL,
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
      } else {
        if (refPlotURL) URL.revokeObjectURL(refPlotURL);
        dispatchProfileComparison({ refErr: true, refPlotURL: '' });
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
    if (fn == 'profileComparisonWithin') {
      dispatchProfileComparison({
        withinSubmitOverlay: true,
        withinErr: '',
        debugR: '',
      });
    } else {
      dispatchProfileComparison({
        refSubmitOverlay: true,
        refErr: '',
        debugR: '',
      });
    }
    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();
        if (fn == 'profileComparisonWithin') {
          dispatchProfileComparison({
            withinSubmitOverlay: false,
            debugR: err,
          });
        } else {
          dispatchProfileComparison({
            refSubmitOverlay: false,
            debugR: err,
          });
        }
      } else {
        const { debugR, output } = await response.json();

        if (fn == 'profileComparisonWithin') {
          dispatchProfileComparison({
            debugR: debugR,
            withinSubmitOverlay: false,
            withinPlotPath: output.plotPath,
          });
          setRPlot(output.plotPath, 'within');
        } else {
          {
            dispatchProfileComparison({
              debugR: debugR,
              refSubmitOverlay: false,
              refPlotPath: output.plotPath,
            });
            setRPlot(output.plotPath, 'refsig');
          }
        }
      }
    } catch (err) {
      dispatchError(err);
      if (fn == 'profileComparisonWithin') {
        dispatchProfileComparison({ withinSubmitOverlay: false });
      } else {
        dispatchProfileComparison({ refSubmitOverlay: false });
      }
    }
  }

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
              <FontAwesomeIcon icon={faPlus} />
            ) : (
              <FontAwesomeIcon icon={faMinus} />
            )}{' '}
            Comparison Within Samples
          </Toggle>
          <Collapse eventKey="0">
            <Body>
              <Form>
                <LoadingOverlay active={withinSubmitOverlay} />
                <div>
                  <Row className="justify-content-center">
                    <Col sm="2">
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
                    <Col sm="4">
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
                    <Col sm="4">
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
                    <Col sm="2" className="m-auto">
                      <Button
                        variant="primary"
                        onClick={() =>
                          calculateR('profileComparisonWithin', {
                            profileType: withinProfileType,
                            sampleName1: withinSampleName1,
                            sampleName2: withinSampleName2,
                          })
                        }
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
                      <div className="d-flex">
                        <a
                          className="px-2 py-1"
                          href={withinPlotURL}
                          download={withinPlotURL.split('/').slice(-1)[0]}
                        >
                          Download Plot
                        </a>
                      </div>
                      <div className="p-2 border rounded">
                        <Row>
                          <Col>
                            <img
                              className="w-100 my-4"
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
              dispatchProfileComparison({
                displayRefSig: !displayRefSig,
              })
            }
          >
            {displayRefSig == true ? (
              <FontAwesomeIcon icon={faPlus} />
            ) : (
              <FontAwesomeIcon icon={faMinus} />
            )}{' '}
            Comparison to Reference Signatures
          </Toggle>
          <Collapse eventKey="1">
            <Body>
              <Form className="my-2">
                <LoadingOverlay active={refSubmitOverlay} />
                <div>
                  <Row className="justify-content-center">
                    <Col sm="2">
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
                    <Col sm="2">
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
                    <Col sm="2">
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
                              className="p-0 font-weight-bold"
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
                    <Col sm="2" className="m-auto">
                      <Button
                        variant="primary"
                        onClick={() =>
                          calculateR('profileComparisonRefSig', {
                            profileType: refProfileType,
                            sampleName: refSampleName,
                            signatureSet: refSignatureSet,
                            compare: refCompare,
                          })
                        }
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
                      <div className="d-flex">
                        <a
                          className="px-2 py-1"
                          href={refPlotURL}
                          download={refPlotURL.split('/').slice(-1)[0]}
                        >
                          Download Plot
                        </a>
                      </div>
                      <div className="p-2 border rounded">
                        <Row>
                          <Col>
                            <img className="w-100 my-4" src={refPlotURL}></img>
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
