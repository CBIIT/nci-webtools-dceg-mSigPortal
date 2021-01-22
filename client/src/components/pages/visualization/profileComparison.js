import React, { useEffect } from 'react';
import {
  Form,
  Row,
  Col,
  Button,
  Popover,
  OverlayTrigger,
  Tab,
  Nav,
} from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faInfoCircle } from '@fortawesome/free-solid-svg-icons';
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
const { Container, Content: TabContent, Pane } = Tab;
const { Item, Link } = Nav;

export default function ProfileComparison({
  submitR,
  getRefSigOptions,
  defaultMatrix,
}) {
  const {
    source,
    study,
    studyOptions,
    cancerType,
    pubExperimentalStrategy,
    pDataOptions,
  } = useSelector((state) => state.visualize);
  const { projectID, matrixList, svgList } = useSelector(
    (state) => state.visualizeResults
  );
  const { profileOptions } = useSelector((state) => state.mutationalProfiles);
  const state = useSelector((state) => state.profileComparison);
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
    filterSignatures,
    refCompare,
    searchFilter,
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
    display,
    withinErr,
    refErr,
    pubErr,
    debugR,
    withinSubmitOverlay,
    refSubmitOverlay,
    pubSubmitOverlay,
  } = state;

  const popover = (
    <Popover id="popover-basic">
      <Title as="h3">{refSignatureSet}</Title>
      <Content>
        <Control
          id="searchSignatures"
          value={searchFilter}
          placeholder="Search Signatures"
          onChange={(e) => {
            const search = e.target.value;
            const signatures = refSignatures.filter(
              (sig) => sig.toLowerCase().indexOf(search.toLowerCase()) > -1
            );
            console.log(signatures);
            dispatchProfileComparison({
              searchFilter: search,
              filterSignatures: signatures,
            });
          }}
        />
        {refSignatures.length > 0 ? (
          <div>
            {filterSignatures.map((signature) => (
              <Button
                key={signature}
                variant="link"
                onClick={() => {
                  console.log('click');
                  let ref = refCompare;
                  if (ref.length) {
                    dispatchProfileComparison({
                      refCompare: (ref += `;1*${signature}`),
                    });
                  } else {
                    dispatchProfileComparison({
                      refCompare: signature,
                    });
                  }
                }}
              >
                {signature}
              </Button>
            ))}
          </div>
        ) : (
          <p>No Signatures Available</p>
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

  useEffect(() => {
    withinPlotPath ? setRPlot(withinPlotPath, 'within') : clearPlot('within');
    refPlotPath ? setRPlot(refPlotPath, 'ref') : clearPlot('ref');
    pubPlotPath ? setRPlot(pubPlotPath, 'pub') : clearPlot('pub');
  }, [withinPlotPath, refPlotPath, pubPlotPath]);

  function setOverlay(type, status) {
    dispatchProfileComparison({ [`${type}SubmitOverlay`]: status });
  }

  async function setRPlot(plotPath, type) {
    try {
      const response = await fetch(`api/results/${projectID}${plotPath}`);
      if (!response.ok) {
        // console.log(await response.json());
      } else {
        const pic = await response.blob();
        const objectURL = URL.createObjectURL(pic);

        if (state[`${type}PlotURL`])
          URL.revokeObjectURL(state[`${type}PlotURL`]);
        dispatchProfileComparison({
          [`${type}PlotURL`]: objectURL,
        });
      }
    } catch (err) {
      dispatchError(err);
    }
  }

  function clearPlot(type) {
    URL.revokeObjectURL(state[`${type}PlotURL`]);
    dispatchProfileComparison({ [`${type}PlotURL`]: '' });
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
        const response = await fetch(`api/getSignaturesR`, {
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
            filterSignatures: signatures,
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

  async function calculateR(type, fn, args) {
    try {
      setOverlay(type, true);

      dispatchProfileComparison({
        [`${type}Err`]: false,
        [`${type}PlotPath`]: '',
        debugR: '',
      });

      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();
        dispatchProfileComparison({ debugR: err });
      } else {
        const { debugR, output, error } = await response.json();

        dispatchProfileComparison({ debugR: debugR });
        if (Object.keys(output).length) {
          dispatchProfileComparison({ [`${type}PlotPath`]: output.plotPath });
        } else {
          dispatchProfileComparison({
            [`${type}Err`]: error || debugR,
            debugR: debugR || error || true,
          });
        }
      }
    } catch (err) {
      dispatchError(err);
    } finally {
      setOverlay(type, false);
    }
  }

  async function getPublicSamples(study, cancerType) {
    try {
      setOverlay('pub', true);
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
    } finally {
      setOverlay('pub', false);
    }
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
    ].sort((a, b) => a - b);

    dispatchProfileComparison({
      userProfileType: profile,
      userMatrixOptions: matrixOptions,
      userMatrixSize: defaultMatrix(profile, matrixOptions),
    });
  }

  let tabs = [
    {
      key: 'within',
      component: (
        <div>
          <Form className="p-3">
            <LoadingOverlay active={withinSubmitOverlay} />
            <Row>
              <Col lg="2">
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
              <Col lg="4">
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
              <Col lg="4">
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
              <Col lg="2" className="d-flex justify-content-end">
                <Button
                  className="mt-auto mb-3"
                  disabled={sampleOptions.length < 2}
                  variant="primary"
                  onClick={() => {
                    if (source == 'user') {
                      calculateR('within', 'profileComparisonWithin', {
                        profileType: withinProfileType,
                        sampleName1: withinSampleName1,
                        sampleName2: withinSampleName2,
                        matrixList: JSON.stringify(
                          matrixList.filter(
                            (matrix) => matrix.Profile_Type == withinProfileType
                          )
                        ),
                      });
                    } else {
                      calculateR('within', 'profileComparisonWithinPublic', {
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
          </Form>
          <div id="pcWithinPlot">
            {withinErr && (
              <div>
                <p>An error has occured. Please verify your input.</p>
                <p>{withinErr}</p>
              </div>
            )}
            {withinPlotURL && (
              <>
                <hr />
                <Plot
                  className="p-3"
                  plotName={withinPlotPath.split('/').slice(-1)[0]}
                  plotURL={withinPlotURL}
                />
              </>
            )}
          </div>
        </div>
      ),
    },
    {
      key: 'reference',
      component: (
        <div>
          <Form className="p-3">
            <LoadingOverlay active={refSubmitOverlay} />
            <Row>
              <Col lg="2">
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
              <Col lg="3">
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
              <Col lg="3">
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
              <Col lg="3">
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
                  />
                  <Text className="text-muted">(Ex. 0.8*SBS5;0.1*SBS1)</Text>
                </Group>
              </Col>
              <Col lg="1" className="d-flex justify-content-end">
                <Button
                  className="mt-auto"
                  style={{ marginBottom: '2.5rem' }}
                  variant="primary"
                  onClick={() => {
                    if (source == 'user') {
                      calculateR('ref', 'profileComparisonRefSig', {
                        profileType: refProfileType,
                        sampleName: refSampleName,
                        signatureSet: refSignatureSet,
                        compare: refCompare,
                        matrixList: JSON.stringify(
                          matrixList.filter(
                            (matrix) => matrix.Profile_Type == withinProfileType
                          )
                        ),
                      });
                    } else {
                      calculateR('ref', 'profileComparisonRefSigPublic', {
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
          </Form>
          <div id="refPlotDownload">
            {refErr && (
              <div>
                <p>An error has occured. Please verify your input.</p>
                <p className="text-danger">{refErr}</p>
              </div>
            )}
            {refPlotURL && (
              <>
                <hr />
                <Plot
                  className="p-3"
                  plotName={refPlotPath.split('/').slice(-1)[0]}
                  plotURL={refPlotURL}
                />
              </>
            )}
          </div>
        </div>
      ),
    },
  ];

  if (source == 'user')
    tabs.push({
      key: 'public',
      component: (
        <div>
          <Form className="p-3">
            <LoadingOverlay active={pubSubmitOverlay} />
            <Row>
              <Col lg="2">
                <Select
                  id="pcUserProfileType"
                  label="Profile Type"
                  value={userProfileType}
                  options={profileOptions}
                  onChange={handleProfile}
                />
              </Col>
              <Col lg="2">
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
              <Col lg="4">
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
              <Col lg="4" />
            </Row>
            <Row>
              <Col lg="2">
                <Select
                  id="pcPubStudy"
                  label="Study"
                  value={pubStudy}
                  options={studyOptions}
                  onChange={handleStudyChange}
                />
              </Col>
              <Col lg="3">
                <Select
                  id="pcPubCancerType"
                  label="Cancer Type"
                  value={pubCancerType}
                  options={pubCancerTypeOptions}
                  onChange={handleCancerChange}
                />
              </Col>
              <Col lg="3">
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
              <Col lg="4" className="d-flex justify-content-end">
                <Button
                  className="mt-auto mb-3"
                  variant="primary"
                  onClick={() =>
                    calculateR('pub', 'profileComparisonPublic', {
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
          </Form>
          <div id="pcPubPlot">
            {pubErr && (
              <div>
                <p>An error has occured. Please verify your input.</p>
              </div>
            )}
            {pubPlotURL && (
              <>
                <hr />
                <Plot
                  className="p-3"
                  plotName={pubPlotPath.split('/').slice(-1)[0]}
                  plotURL={pubPlotURL}
                />
              </>
            )}
          </div>
        </div>
      ),
    });
  return (
    <div>
      <Container
        transition={false}
        className="mt-2"
        defaultActiveKey={display}
        activeKey={display}
        onSelect={(tab) => dispatchProfileComparison({ display: tab })}
      >
        <Nav variant="tabs">
          <Item>
            <Link eventKey="within" as="button" className="outline-none">
              <strong>PC Within Samples</strong>
            </Link>
          </Item>
          <Item>
            <Link eventKey="reference" as="button" className="outline-none">
              <strong>PC to Reference Signatures</strong>
            </Link>
          </Item>
          {source == 'user' && (
            <Item>
              <Link eventKey="public" as="button" className="outline-none">
                <strong>PC to Public Data</strong>
              </Link>
            </Item>
          )}
        </Nav>
        <TabContent
          className={`bg-white tab-pane-bordered rounded-0 d-block`}
          style={{ overflowX: 'auto' }}
        >
          {tabs.map(({ key, component }) => (
            <Pane key={key} eventKey={key} className="border-0">
              {component}
            </Pane>
          ))}
        </TabContent>
      </Container>
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
