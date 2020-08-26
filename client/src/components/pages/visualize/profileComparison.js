import React, { useEffect } from 'react';
import {
  Form,
  Row,
  Col,
  Button,
  Popover,
  OverlayTrigger,
} from 'react-bootstrap';
import { Typeahead } from 'react-bootstrap-typeahead';
import { useSelector } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faInfoCircle } from '@fortawesome/free-solid-svg-icons';
import {
  dispatchError,
  dispatchProfileComparison,
} from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import 'react-bootstrap-typeahead/css/Typeahead.css';

const { Group, Label, Control, Text } = Form;
const { Title, Content } = Popover;

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
    debugR,
    displayDebug,
    withinSubmitOverlay,
    refSubmitOverlay,
  } = useSelector((state) => state.profileComparison);

  // load r plots after they are recieved
  // useEffect(() => {
  //   if (displayTab == 'profileComparison') {
  //     if (withinPlotPath && !withinPlotURL && !withinSubmitOverlay)
  //       setRPlot(withinPlotPath, 'within');
  //     if (refPlotPath && !refPlotURL && !refSubmitOverlay)
  //       setRPlot(refPlotPath, 'refsig');
  //   }
  // }, [withinPlotPath, refPlotPath, displayTab]);

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
        setOverlay(type, false);
      } else {
        const pic = await response.blob();
        const objectURL = URL.createObjectURL(pic);

        if (type == 'within') {
          if (withinPlotURL.length) URL.revokeObjectURL(withinPlotURL);
          dispatchProfileComparison({
            withinPlotURL: objectURL,
          });
        } else {
          if (refPlotURL.length) URL.revokeObjectURL(refPlotURL);
          dispatchProfileComparison({
            refPlotURL: objectURL,
          });
        }
        setOverlay(type, false);
      }
    } catch (err) {
      dispatchError(err);
      setOverlay(type, false);
    }
  }

  // get Signature Reference Sets for dropdown options
  async function getSignatureSet(profileType) {
    if (profileType && profileType.length) {
      dispatchProfileComparison({ refSubmitOverlay: true });
      try {
        const response = await getRefSigOptions(profileType);

        if (response.ok) {
          const signatureSetOptions = await response.json();

          dispatchProfileComparison({
            refProfileType: profileType,
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
    if (signatureSetName && signatureSetName.length) {
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
        debugR: '',
      });
    } else {
      dispatchProfileComparison({
        refSubmitOverlay: true,
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
      <Form>
        <LoadingOverlay active={withinSubmitOverlay} />
        <Label>
          <Button
            variant="link"
            className="p-0 font-weight-bold"
            onClick={() =>
              dispatchProfileComparison({
                displayWithin: !displayWithin,
              })
            }
          >
            Comparison Within Samples
          </Button>
        </Label>
        <div
          className="border rounded p-2"
          style={{ display: displayWithin ? 'block' : 'none' }}
        >
          <Row className="justify-content-center">
            <Col sm="2">
              <Group controlId="profileTypeWithin">
                <Label>Profile Type</Label>
                <Control
                  as="select"
                  value={withinProfileType}
                  onChange={(e) =>
                    dispatchProfileComparison({
                      withinProfileType: e.target.value,
                    })
                  }
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
            <Col sm="4">
              <Label>Sample Name 1</Label>
              <Typeahead
                clearButton
                id="selectSampleName1"
                placeholder="Select Sample Name 1"
                selected={[withinSampleName1]}
                options={nameOptions}
                onChange={(name) =>
                  dispatchProfileComparison({
                    withinSampleName1: name[0] || '',
                  })
                }
              />
              {/* <Control
                as="select"
                value={withinSampleName1}
                onChange={(e) =>
                  dispatchProfileComparison({
                    withinSampleName1: e.target.value,
                  })
                }
                custom
              >
                <option value="0" disabled>
                  Select
                </option>
                {nameOptions.map((name, index) => {
                  return (
                    <option key={index} value={name}>
                      {name}
                    </option>
                  );
                })}
              </Control> */}
            </Col>
            <Col sm="4">
              <Label>Sample Name 2</Label>
              <Typeahead
                clearButton
                id="selectSampleName2"
                placeholder="Select Sample Name 2"
                selected={[withinSampleName2]}
                options={nameOptions}
                onChange={(name) =>
                  dispatchProfileComparison({
                    withinSampleName2: name[0] || '',
                  })
                }
              />
              {/* <Control
                as="select"
                value={withinSampleName2}
                onChange={(e) =>
                  dispatchProfileComparison({
                    withinSampleName2: e.target.value,
                  })
                }
                custom
              >
                <option value="0" disabled>
                  Select
                </option>
                {nameOptions.map((name, index) => {
                  return (
                    <option key={index} value={name}>
                      {name}
                    </option>
                  );
                })}
              </Control> */}
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

          <div
            id="pcWithinPlot"
            style={{ display: withinPlotURL.length ? 'block' : 'none' }}
          >
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
                  <img className="w-100 my-4" src={withinPlotURL}></img>
                </Col>
              </Row>
            </div>
          </div>
        </div>
      </Form>

      <Form className="my-2">
        <LoadingOverlay active={refSubmitOverlay} />
        <Label>
          <Button
            variant="link"
            className="p-0 font-weight-bold"
            onClick={() =>
              dispatchProfileComparison({
                displayRefSig: !displayRefSig,
              })
            }
          >
            Comparison to Reference Signatures
          </Button>
        </Label>
        <div
          className="border rounded p-2"
          style={{ display: displayRefSig ? 'block' : 'none' }}
        >
          <Row className="justify-content-center">
            <Col sm="2">
              <Group controlId="profileTypeRefSig">
                <Label>Profile Type</Label>
                <Control
                  as="select"
                  value={refProfileType}
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
            <Col sm="2">
              <Group controlId="sampleNameRefSig">
                <Label>Sample Name</Label>
                <Typeahead
                  clearButton
                  id="selectSampleNameRef"
                  placeholder="Select Sample Name"
                  selected={[refSampleName]}
                  options={nameOptions}
                  onChange={(name) =>
                    dispatchProfileComparison({
                      refSampleName: name[0] || '',
                    })
                  }
                />
                {/* <Control
                  as="select"
                  value={refSampleName}
                  onChange={(e) => {
                    dispatchProfileComparison({
                      refSampleName: e.target.value,
                    });
                  }}
                  custom
                > 
                  {nameOptions.map((name, index) => {
                    return (
                      <option key={index} value={name}>
                        {name}
                      </option>
                    );
                  })}
                </Control>*/}
              </Group>
            </Col>
            <Col sm="4">
              <Group controlId="signatureSet">
                <Label>Reference Signature Set</Label>
                <Control
                  disabled={!refSignatureSetOptions.length}
                  as="select"
                  value={refSignatureSet}
                  onChange={(e) => {
                    dispatchProfileComparison({
                      refSignatureSet: e.target.value,
                    });
                    getSignatures(refProfileType, e.target.value);
                  }}
                  custom
                >
                  <option value="0">Select</option>
                  {refSignatureSetOptions.map((signatureSet, index) => {
                    return (
                      <option key={index} value={signatureSet}>
                        {signatureSet}
                      </option>
                    );
                  })}
                </Control>
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
                    <Button variant="link" className="p-0 font-weight-bold">
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
                <Text className="text-muted">(Ex. 0.8*SBS5;0.1*SBS1)</Text>
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

          <div
            id="refPlotDownload"
            style={{ display: refPlotURL.length ? 'block' : 'none' }}
          >
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
      </Form>

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
