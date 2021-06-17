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
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faInfoCircle } from '@fortawesome/free-solid-svg-icons';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';
import Debug from '../../controls/debug/debug';
import Select from '../../controls/select/select';
import {
  value2d,
  filter2d,
  unique2d,
  defaultMatrix,
} from '../../../services/utils';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../services/store/visualization';
import { actions as modalActions } from '../../../services/store/modal';

const actions = { ...visualizationActions, ...modalActions };
const { Group, Label, Control, Text } = Form;
const { Title, Content } = Popover;
const { Container, Content: TabContent, Pane } = Tab;
const { Item, Link } = Nav;

export default function ProfileComparison({ submitR, getRefSigOptions }) {
  const dispatch = useDispatch();
  const visualization = useSelector((state) => state.visualization);
  const mergeProfileComparison = (state) =>
    dispatch(actions.mergeVisualization({ profileComparison: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const {
    source,
    study,
    studyOptions,
    cancerType,
    pubExperimentalStrategy,
    pDataOptions,
  } = visualization.visualize;
  const { projectID, matrixList, svgList } = visualization.results;
  const { profileOptions } = visualization.mutationalProfiles;
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
  } = visualization.profileComparison;

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

            mergeProfileComparison({
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
                  let ref = refCompare;
                  if (ref.length) {
                    mergeProfileComparison({
                      refCompare: (ref += `;1*${signature}`),
                    });
                  } else {
                    mergeProfileComparison({
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
    mergeProfileComparison({ [`${type}SubmitOverlay`]: status });
  }

  async function setRPlot(plotPath, type) {
    try {
      const response = await fetch(`api/results/${projectID}${plotPath}`);
      if (!response.ok) {
        // console.log(await response.json());
      } else {
        const pic = await response.blob();
        const objectURL = URL.createObjectURL(pic);

        if (visualization.profileComparison[`${type}PlotURL`])
          URL.revokeObjectURL(
            visualization.profileComparison[`${type}PlotURL`]
          );
        mergeProfileComparison({
          [`${type}PlotURL`]: objectURL,
        });
      }
    } catch (err) {
      mergeError(err.message);
    }
  }

  function clearPlot(type) {
    URL.revokeObjectURL(visualization.profileComparison[`${type}PlotURL`]);
    mergeProfileComparison({ [`${type}PlotURL`]: '' });
  }

  // get Signature Reference Sets for dropdown options
  async function getSignatureSet(profileType) {
    if (profileType) {
      mergeProfileComparison({ refSubmitOverlay: true });
      try {
        const response = await getRefSigOptions(profileType);

        if (response.ok) {
          const signatureSetOptions = await response.json();

          mergeProfileComparison({
            refSignatureSetOptions: signatureSetOptions,
            refSignatureSet: signatureSetOptions[0],
            refSubmitOverlay: false,
          });
          getSignatures(profileType, signatureSetOptions[0]);
        } else {
          mergeError(await response.json());
          mergeProfileComparison({ refSubmitOverlay: false });
        }
      } catch (err) {
        mergeError(err.message);
        mergeProfileComparison({ refSubmitOverlay: false });
      }
    }
  }

  // get signature options for compare
  async function getSignatures(profileType, signatureSetName) {
    if (signatureSetName) {
      mergeProfileComparison({ refSubmitOverlay: true });
      try {
        const response = await fetch(`api/visualizationData`, {
          method: 'POST',
          headers: {
            Accept: 'application/json',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            fn: 'getSignatureNames',
            args: {
              profileType: profileType,
              signatureSetName: signatureSetName,
            },
          }),
        });

        if (response.ok) {
          const signatures = await response.json();

          mergeProfileComparison({
            refSignatures: signatures,
            filterSignatures: signatures,
            refCompare: signatures[0],
            refSubmitOverlay: false,
          });
        } else {
          mergeError(await response.json());
          mergeProfileComparison({ refSubmitOverlay: false });
        }
      } catch (err) {
        mergeError(err.message);
        mergeProfileComparison({ refSubmitOverlay: false });
      }
    }
  }

  async function calculateR(type, fn, args) {
    try {
      setOverlay(type, true);

      mergeProfileComparison({
        [`${type}Err`]: false,
        [`${type}PlotPath`]: '',
        debugR: '',
      });

      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();
        mergeProfileComparison({ debugR: err });
      } else {
        const { debugR, output, error } = await response.json();

        mergeProfileComparison({ debugR: debugR });
        if (Object.keys(output).length) {
          mergeProfileComparison({ [`${type}PlotPath`]: output.plotPath });
        } else {
          mergeProfileComparison({
            [`${type}Err`]: error || debugR,
            debugR: debugR || error || true,
          });
        }
      }
    } catch (err) {
      mergeError(err.message);
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
        const pubSamples = unique2d('Sample', svgList.columns, svgList.data);

        mergeProfileComparison({
          pubSampleOptions: pubSamples,
          pubSampleName: pubSamples[0],
        });
      } else {
        mergeProfileComparison({
          pubSampleOptions: ['NA'],
          pubErr: 'Error retrieving pub data samples',
        });
      }
    } catch (err) {
      mergeError(err.message);
      mergeProfileComparison({
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

    mergeProfileComparison({
      pubStudy: study,
      pubCancerType: cancerTypeOptions[0],
      pubCancerTypeOptions: cancerTypeOptions,
    });
    getPublicSamples(study, cancerTypeOptions[0]);
  }

  function handleCancerChange(cancer) {
    mergeProfileComparison({
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

    mergeProfileComparison({
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
          <div className="p-3">
            <p>
              Input a “Profile Type” and two “Sample Names” to generate the
              mutational profile of each sample, as well as the difference
              between the two mutational profiles.
            </p>
          </div>
          <hr />
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
                    mergeProfileComparison({
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
                    mergeProfileComparison({
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
                    mergeProfileComparison({
                      withinSampleName2: name,
                    });
                  }}
                />
              </Col>
              <Col lg="2" className="d-flex">
                <Button
                  className="ml-auto mb-auto"
                  disabled={sampleOptions.length < 2}
                  variant="primary"
                  onClick={() => {
                    if (source == 'user') {
                      calculateR('within', 'profileComparisonWithin', {
                        profileType: withinProfileType,
                        sampleName1: withinSampleName1,
                        sampleName2: withinSampleName2,
                        matrixFile: value2d(
                          filter2d(
                            [
                              withinProfileType,
                              defaultMatrix(withinProfileType, [
                                '96',
                                '78',
                                '83',
                              ]),
                            ],
                            matrixList.data
                          )[0],
                          'Path',
                          matrixList.columns
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
                  downloadName={withinPlotPath.split('/').slice(-1)[0]}
                  plotURL={withinPlotURL}
                />
                <div className="p-3">
                  <p>
                    The plot generated displays the mutational profile for the
                    two samples selected, and the difference between them. Also
                    at the top of the plot are measurements for RSS and cosine
                    similarity.
                  </p>
                  <p>
                    RSS is the Residual Sum of Squares. It measures the
                    discrepancy between two profiles. Cosine similarity is how
                    similar the mutational profiles are to one another. For
                    additional information about RSS and cosine similarity,
                    click <a href="#faq">here.</a>
                  </p>
                </div>
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
          <div className="p-3">
            <p>
              Input a “Profile Type”, a “Sample Name”, a “Reference Signature
              Set”, and a “Compare Signatures” parameter to generate the
              mutational profile of the sample, the signature from the reference
              set, and the difference between the two mutational profiles.
            </p>
          </div>
          <hr />
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
                    mergeProfileComparison({
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
                    mergeProfileComparison({
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
                    mergeProfileComparison({
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
                      mergeProfileComparison({
                        refCompare: e.target.value,
                      });
                    }}
                  />
                  <Text className="text-muted">(Ex. 0.8*SBS5;0.1*SBS1)</Text>
                </Group>
              </Col>
              <Col lg="1" className="d-flex">
                <Button
                  className="ml-auto mb-auto"
                  style={{ marginBottom: '2.5rem' }}
                  variant="primary"
                  onClick={() => {
                    if (source == 'user') {
                      calculateR('ref', 'profileComparisonRefSig', {
                        profileType: refProfileType,
                        sampleName: refSampleName,
                        signatureSet: refSignatureSet,
                        compare: refCompare,
                        matrixFile: value2d(
                          filter2d(
                            [
                              withinProfileType,
                              defaultMatrix(withinProfileType, [
                                '96',
                                '78',
                                '83',
                              ]),
                            ],
                            matrixList.data
                          )[0],
                          'Path',
                          matrixList.columns
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
                  downloadName={refPlotPath.split('/').slice(-1)[0]}
                  plotURL={refPlotURL}
                />
                <div className="p-3">
                  <p>
                    The plot generated displays the mutational profile for the
                    sample, the signature selected from the reference signature
                    set, and the difference between them. Also at the top of the
                    plot are measurements for RSS and cosine similarity.
                  </p>
                  <p>
                    RSS is the Residual Sum of Squares. It measures the
                    discrepancy between two profiles. Cosine similarity is how
                    similar the mutational profiles are to one another. For
                    additional information about RSS and cosine similarity click{' '}
                    <a href="#faq">here.</a>
                  </p>
                </div>
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
          <div className="p-3">
            <p>
              Input the “Profile Type”, “Matrix Size”, “Sample Name” (from your
              input data), the public “Study”, “Cancer Type”, and “Public Sample
              Name”.
            </p>
          </div>
          <hr />
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
                    mergeProfileComparison({
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
                    mergeProfileComparison({
                      userSampleName: name,
                    });
                  }}
                />
              </Col>
              <Col />
              <Col lg="2" className="d-flex">
                <Button
                  className="ml-auto mb-auto"
                  variant="primary"
                  onClick={() =>
                    calculateR('pub', 'profileComparisonPublic', {
                      profileName: userProfileType + userMatrixSize,
                      matrixFile: value2d(
                        filter2d(
                          [userProfileType, userMatrixSize],
                          matrixList.data
                        )[0],
                        'Path',
                        matrixList.columns
                      ),
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
                    mergeProfileComparison({
                      pubSampleName: name,
                    });
                  }}
                />
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
                  downloadName={pubPlotPath.split('/').slice(-1)[0]}
                  plotURL={pubPlotURL}
                />
                <div className="p-3">
                  <p>
                    The plot generated displays the mutational profile of the
                    sample from the input data, the sample from the public data
                    selected, and the difference between them. Also at the top
                    of the plot are measurements for RSS and cosine similarity.
                  </p>
                  <p>
                    RSS is the Residual Sum of Squares. It measures the
                    discrepancy between two profiles. Cosine similarity is how
                    similar the mutational profiles are to one another. For
                    additional information about RSS and cosine similarity,
                    click <a href="#faq">here.</a>
                  </p>
                </div>
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
        onSelect={(tab) => mergeProfileComparison({ display: tab })}
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
          <div className="px-3 pt-3 pb-0">
            <p className="m-0">
              Below you can observe mutational profile comparisons between
              samples (PC Within Samples), between a sample and a signature from
              a reference signature set (PC to Reference Signatures), or between
              user data and public data (PC to Public Data).
            </p>
          </div>
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
