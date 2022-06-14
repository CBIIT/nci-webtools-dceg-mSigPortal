import React, { useEffect, useState } from 'react';
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
import CustomSelect from '../../controls/select/select-old';
import Description from '../../controls/description/description';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../services/store/visualization';
import { actions as modalActions } from '../../../services/store/modal';
import { defaultMatrix } from '../../../services/utils';
import { NavHashLink } from 'react-router-hash-link';

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
    projectID,
    matrixList,
    svgList,
    profileOptions,
  } = visualization.state;

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
    withinTxtPath,
    refTxtPath,
    pubTxtPath,
    display,
    withinErr,
    refErr,
    pubErr,
    debugR,
    withinSubmitOverlay,
    refSubmitOverlay,
    pubSubmitOverlay,
  } = visualization.profileComparison;

  const [invalidSignature, setInvalid] = useState(false);

  const popover = (
    <Popover id="popover-basic" style={{ minWidth: '400px' }}>
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
    if (refProfileType && refSignatureSet && !refSignatures.length)
      getSignatures(refProfileType, refSignatureSet);
  }, [refProfileType, refSignatureSet]);

  // set initial study
  useEffect(() => {
    if (!pubStudy && study) handleStudyChange(study);
  }, []);

  // get public data samples
  useEffect(() => {
    if (!pubSampleOptions.length && source == 'user' && pubStudy)
      getPublicSamples(pubStudy, pubCancerType);
  }, [pDataOptions, pubStudy]);

  function setOverlay(type, status) {
    mergeProfileComparison({ [`${type}SubmitOverlay`]: status });
  }

  // get Signature Reference Sets for dropdown options
  async function getSignatureSet(profileType) {
    if (profileType) {
      mergeProfileComparison({ refSubmitOverlay: true });
      try {
        const response = await getRefSigOptions(profileType);

        if (response.ok) {
          const { output: signatureSetOptions } = await response.json();

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
        const response = await fetch(`web/visualizationWrapper`, {
          method: 'POST',
          headers: {
            Accept: 'application/json',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            fn: 'getSignaturesR',
            args: {
              profileType: profileType,
              signatureSetName: signatureSetName,
            },
          }),
        });

        if (response.ok) {
          const { output: signatures } = await response.json();

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
        [`${type}TxtPath`]: '',
        debugR: '',
      });

      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();
        mergeProfileComparison({ debugR: err });
      } else {
        const { output } = await response.json();

        mergeProfileComparison({ debugR: debugR });
        if (output.plotPath) {
          mergeProfileComparison({
            [`${type}PlotPath`]: output.plotPath,
            [`${type}TxtPath`]: output.txtPath,
          });
        } else {
          mergeProfileComparison({
            [`${type}Err`]: output.error || output.uncaughtError || true,
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

      const response = await fetch(`web/visualizationWrapper`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ fn: 'getPublicData', args }),
      });

      if (response.ok) {
        const { output: svgList } = await response.json();
        const pubSamples = [...new Set(svgList.map(({ Sample }) => Sample))];

        mergeProfileComparison({
          pubSampleOptions: pubSamples,
          pubSampleName: pubSamples[0],
          pubErr: false,
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
          .filter((row) => row.Profile_Type == profile)
          .map(({ Matrix_Size }) => Matrix_Size)
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
          <p className="p-3 m-0">
            Input a [Profile Type] and two sample names ([Sample Name 1],
            [Sample Name 2]) to generate the mutational profile of each sample,
            as well as the difference between the two mutational profiles.
          </p>

          <hr />
          <Form className="p-3">
            <LoadingOverlay active={withinSubmitOverlay} />
            <Row>
              <Col lg="auto">
                <CustomSelect
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
              <Col lg="auto">
                <CustomSelect
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
              <Col lg="auto">
                <CustomSelect
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
              <Col lg="auto" className="d-flex">
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
                        matrixFile: matrixList.filter(
                          (row) =>
                            row.Profile_Type == withinProfileType &&
                            row.Matrix_Size ==
                              defaultMatrix(withinProfileType, [
                                '96',
                                '78',
                                '83',
                              ])
                        )[0].Path,
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
              <>
                <hr />
                <div className="p-3">
                  <p>An error has occured. Please verify your input.</p>
                  <p>{withinErr}</p>
                </div>
              </>
            )}
            {withinPlotPath && (
              <>
                <hr />
                <Plot
                  className="p-3"
                  downloadName={withinPlotPath.split('/').slice(-1)[0]}
                  plotPath={'web/results/' + withinPlotPath}
                  txtPath={`web/results/${withinTxtPath}`}
                />
                <div className="p-3">
                  <p>
                    The plot above shows the mutational profiles of two selected
                    samples, as well as the difference between them. The text at
                    the top of the plot indicates the profile similarity
                    calculated using Residual Sum of Squares (RSS) and cosine
                    similarity methods.
                  </p>
                  <p>
                    RSS measures the discrepancy between two mutational
                    profiles. Cosine similarity measures how similar two
                    mutational profiles are. For example, two identical
                    mutational profiles will have RSS = 0 and Cosine similarity
                    = 1. For additional information about RSS and cosine
                    similarity, click{' '}
                    <NavHashLink to="/faq#cosine-similarity">here</NavHashLink>.
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
          <p className="p-3 m-0">
            Input a [Profile Type], [Sample Name], [Reference Signature Set],
            and [Compare Signatures] parameter to generate the mutational
            profile of the sample, the signature profile from the reference set,
            and the difference between two mutational profiles. The [Compare
            Signatures] can be a single reference signature name or a combined
            multiple reference signature name with different contributions (see
            example below [Compare Signatures]).
          </p>
          <hr />
          <Form className="p-3">
            <LoadingOverlay active={refSubmitOverlay} />
            <Row>
              <Col lg="auto">
                <CustomSelect
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
              <Col lg="auto">
                <CustomSelect
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
              <Col lg="auto">
                <CustomSelect
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
              <Col lg="auto">
                <Group controlId="refCompare">
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
                    isInvalid={invalidSignature}
                    onChange={(e) => {
                      mergeProfileComparison({
                        refCompare: e.target.value,
                      });
                    }}
                  />
                  <Text className="text-muted">(Ex. 0.8*SBS5;0.1*SBS1)</Text>
                  <Form.Control.Feedback type="invalid">
                    Enter a valid signature. Click info icon for options.
                  </Form.Control.Feedback>
                </Group>
              </Col>
              <Col lg="auto" className="d-flex">
                <Button
                  className="mt-auto"
                  style={{ marginBottom: '2.5rem' }}
                  variant="primary"
                  onClick={() => {
                    if (refCompare) {
                      setInvalid(false);
                      if (source == 'user') {
                        calculateR('ref', 'profileComparisonRefSig', {
                          profileType: refProfileType,
                          sampleName: refSampleName,
                          signatureSet: refSignatureSet,
                          compare: refCompare,
                          matrixFile: matrixList.filter(
                            (row) =>
                              row.Profile_Type == withinProfileType &&
                              row.Matrix_Size ==
                                defaultMatrix(withinProfileType, [
                                  '96',
                                  '78',
                                  '83',
                                ])
                          )[0].Path,
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
                    } else {
                      setInvalid(true);
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
              <>
                <hr />
                <div className="p-3">
                  <p>An error has occured. Please verify your input.</p>
                  <p>{refErr}</p>
                </div>
              </>
            )}
            {refPlotPath && (
              <>
                <hr />
                <Plot
                  className="p-3"
                  downloadName={refPlotPath.split('/').slice(-1)[0]}
                  plotPath={'web/results/' + refPlotPath}
                  txtPath={`web/results/${refTxtPath}`}
                />
                <div className="p-3">
                  <p>
                    The plot above shows the mutational profiles of a selected
                    sample, the signature from the selected reference signature
                    set, and the difference between them. The text at the top of
                    the plot indicates the profile similarity calculated using
                    Residual Sum of Squares (RSS) and cosine similarity methods.
                  </p>
                  <p>
                    RSS measures the discrepancy between two mutational
                    profiles. Cosine similarity measures how similar two
                    mutational profiles are. For example, two identical
                    mutational profiles will have RSS = 0 and Cosine similarity
                    = 1. For additional information about RSS and cosine
                    similarity, click{' '}
                    <NavHashLink to="/faq#cosine-similarity">here</NavHashLink>.
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
          <p className="p-3 m-0">
            Input a [Profile Type], [Matrix Size], [Sample Name] (from your
            input data), [Study], [Cancer Type], and [Public Sample Name] to
            generate the mutational profile of the input sample, the sample from
            the selected public data, and the difference between them.
          </p>
          <hr />
          <Form className="p-3">
            <LoadingOverlay active={pubSubmitOverlay} />
            <Row>
              <Col lg="auto">
                <CustomSelect
                  id="pcUserProfileType"
                  label="Profile Type"
                  value={userProfileType}
                  options={profileOptions}
                  onChange={(profile) => handleProfile(profile)}
                />
              </Col>
              <Col lg="auto">
                <CustomSelect
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
              <Col lg="auto">
                <CustomSelect
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
              <Col lg="auto">
                <CustomSelect
                  id="pcPubStudy"
                  label="Study"
                  value={pubStudy}
                  options={studyOptions}
                  onChange={handleStudyChange}
                />
              </Col>
              <Col lg="auto">
                <CustomSelect
                  id="pcPubCancerType"
                  label="Cancer Type or Group"
                  value={pubCancerType}
                  options={pubCancerTypeOptions}
                  onChange={handleCancerChange}
                />
              </Col>
              <Col lg="auto">
                <CustomSelect
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
              <Col lg="auto" className="d-flex">
                <Button
                  className="mt-auto mb-3"
                  variant="primary"
                  onClick={() =>
                    calculateR('pub', 'profileComparisonPublic', {
                      profileName: userProfileType + userMatrixSize,
                      matrixFile: matrixList.filter(
                        (row) =>
                          row.Profile_Type == userProfileType &&
                          row.Matrix_Size == userMatrixSize
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
              <>
                <hr />
                <div className="p-3">
                  <p>{pubErr}</p>
                </div>
              </>
            )}
            {pubPlotPath && (
              <>
                <hr />
                <Plot
                  className="p-3"
                  downloadName={pubPlotPath.split('/').slice(-1)[0]}
                  plotPath={'web/results/' + pubPlotPath}
                  txtPath={`web/results/${pubTxtPath}`}
                />
                <div className="p-3">
                  <p>
                    The plot above shows the mutational profile of the input
                    sample, the sample from the selected public data, and the
                    difference between them. The text on the top of the plot
                    indicates the profile similarity calculated using Residual
                    Sum of Squares (RSS) and cosine similarity methods.
                  </p>
                  <p>
                    RSS measures the discrepancy between two mutational
                    profiles. Cosine similarity measures how similar two
                    mutational profiles are. For example, two identical
                    mutational profiles will have RSS = 0 and Cosine similarity
                    = 1. For additional information about RSS and cosine
                    similarity, click{' '}
                    <NavHashLink to="/faq#cosine-similarity">here</NavHashLink>.
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
      <p className="bg-white border rounded p-3 mb-1">
        Below you can perform mutational profile comparison analyses between
        samples (PC Between Samples), between a sample and a signature from a
        selected reference signature set (PC to Reference Signatures), or
        between user data and public data (PC to Public Data).
      </p>
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
              <strong>PC Between Samples</strong>
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
