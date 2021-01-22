import React, { useState, useEffect } from 'react';
import { Form, Row, Col, Button, Tab, Nav } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchCosineSimilarity,
} from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';
import Debug from '../../controls/debug/debug';
import Select from '../../controls/select/select';

const { Container, Content, Pane } = Tab;
const { Item, Link } = Nav;

export default function CosineSimilarity({
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
  const { profileOptions } = useSelector((state) => state.mutationalProfiles);
  const { projectID, matrixList, svgList } = useSelector(
    (state) => state.visualizeResults
  );
  const state = useSelector((state) => state.cosineSimilarity);
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
    display,
    withinErr,
    refErr,
    pubErr,
    withinSubmitOverlay,
    refSubmitOverlay,
    pubSubmitOverlay,
    debugR,
  } = state;

  // check for multiple sample input and disable params if true
  const [multiSample, setMultiSample] = useState(false);
  useEffect(() => {
    if (svgList.length) {
      if (source == 'user') {
        const samples = [
          ...new Set(
            svgList.map((plot) => {
              if (plot.Filter != 'NA')
                return `${plot.Sample_Name}@${plot.Filter}`;
              else return plot.Sample_Name;
            })
          ),
        ];

        if (samples.length > 1) setMultiSample(true);
        else setMultiSample(false);
      } else {
        setMultiSample(true);
      }
    }
  }, [svgList]);

  useEffect(() => {
    withinPlotPath ? setRPlot(withinPlotPath, 'within') : clearPlot('within');
    refPlotPath ? setRPlot(refPlotPath, 'ref') : clearPlot('ref');
    pubPlotPath ? setRPlot(pubPlotPath, 'pub') : clearPlot('pub');
  }, [withinPlotPath, refPlotPath, pubPlotPath]);

  function setOverlay(type, display) {
    dispatchCosineSimilarity({ [`${type}SubmitOverlay`]: display });
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
        dispatchCosineSimilarity({
          [`${type}PlotURL`]: objectURL,
        });
      }
    } catch (err) {
      dispatchError(err);
    }
  }

  function clearPlot(type) {
    URL.revokeObjectURL(state[`${type}PlotURL`]);
    dispatchCosineSimilarity({ [`${type}PlotURL`]: '' });
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

  function setOverlay(type, status) {
    dispatchCosineSimilarity({ [`${type}SubmitOverlay`]: status });
  }

  async function calculateR(type, fn, args) {
    try {
      setOverlay(type, true);

      dispatchCosineSimilarity({
        [`${type}Err`]: false,
        [`${type}PlotPath`]: '',
        [`${type}TxtPath`]: '',
        debugR: '',
      });

      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();
        dispatchCosineSimilarity({ debugR: err });
      } else {
        const { debugR, output, error } = await response.json();

        if (Object.keys(output).length) {
          dispatchCosineSimilarity({
            [`${type}PlotPath`]: output.plotPath,
            [`${type}TxtPath`]: output.txtPath,
          });
        } else {
          dispatchCosineSimilarity({
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

  function handleWithinProfileType(profileType) {
    if (source == 'user') {
      const withinMatrixOptions = [
        ...new Set(
          matrixList
            .filter((matrix) => matrix.Profile_Type == profileType)
            .map((matrix) => matrix.Matrix_Size)
        ),
      ].sort((a, b) => a - b);

      dispatchCosineSimilarity({
        withinProfileType: profileType,
        withinMatrixSize: defaultMatrix(profileType, withinMatrixOptions),
        withinMatrixOptions: withinMatrixOptions,
      });
    } else {
      const withinMatrixOptions = [
        ...new Set(
          svgList
            .filter((plot) => plot.Profile.indexOf(profileType) > -1)
            .map((plot) => plot.Profile.match(/\d+/gi)[0])
        ),
      ].sort((a, b) => a - b);
      dispatchCosineSimilarity({
        withinProfileType: profileType,
        withinMatrixSize: defaultMatrix(profileType, withinMatrixOptions),
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
    ].sort((a, b) => a - b);

    dispatchCosineSimilarity({
      userProfileType: profileType,
      userMatrixSize: defaultMatrix(profileType, userMatrixOptions),
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
                  disabled={!multiSample}
                  id="csProfileType"
                  label="Profile Type"
                  value={withinProfileType}
                  options={profileOptions}
                  onChange={handleWithinProfileType}
                />
              </Col>
              <Col lg="2">
                <Select
                  disabled={!multiSample}
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
              <Col lg="6" />
              <Col lg="2" className="d-flex justify-content-end">
                <Button
                  className="mt-auto mb-3"
                  disabled={!multiSample}
                  variant="primary"
                  onClick={() => {
                    if (source == 'user') {
                      calculateR('within', 'cosineSimilarityWithin', {
                        matrixFile: matrixList.filter(
                          (path) =>
                            path.Profile_Type == withinProfileType &&
                            path.Matrix_Size == withinMatrixSize
                        )[0].Path,
                      });
                    } else {
                      calculateR('within', 'cosineSimilarityWithinPublic', {
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
            {!multiSample && (
              <Row>
                <Col>Unavailable - More than one Sample Required</Col>
              </Row>
            )}
          </Form>

          <div id="withinPlot">
            {withinErr && (
              <div className="p-3">
                <p>An error has occured. Please verify your input.</p>
              </div>
            )}
            {withinPlotURL && (
              <>
                <hr />
                <Plot
                  className="p-3"
                  plotName={withinPlotPath.split('/').slice(-1)[0]}
                  plotURL={withinPlotURL}
                  txtPath={projectID + withinTxtPath}
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
            <Row className="justify-content-center">
              <Col lg="2">
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
              <Col lg="3">
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
              <Col lg="5" />
              <Col lg="2" className="d-flex justify-content-end">
                <Button
                  className="mt-auto mb-3"
                  variant="primary"
                  onClick={() => {
                    if (source == 'user') {
                      calculateR('ref', 'cosineSimilarityRefSig', {
                        profileType: refProfileType,
                        signatureSet: refSignatureSet,
                        matrixList: JSON.stringify(
                          matrixList.filter(
                            (matrix) => matrix.Profile_Type == refProfileType
                          )
                        ),
                      });
                    } else {
                      calculateR('ref', 'cosineSimilarityRefSigPublic', {
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
          </Form>
          <div id="refPlot">
            {refErr && (
              <div className="p-3">
                <p>An error has occured. Please verify your input.</p>
              </div>
            )}
            {refPlotURL && (
              <>
                <hr />
                <Plot
                  className="p-3"
                  plotName={refPlotPath.split('/').slice(-1)[0]}
                  plotURL={refPlotURL}
                  txtPath={projectID + refTxtPath}
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
            <Row className="justify-content-center">
              <Col lg="2">
                <Select
                  id="csUserProfileType"
                  label="Profile Type"
                  value={userProfileType}
                  options={profileOptions}
                  onChange={handlePublicProfileType}
                />
              </Col>
              <Col lg="2">
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
              <Col lg="2">
                <Select
                  id="csPubStudy"
                  label="Study"
                  value={pubStudy}
                  options={studyOptions}
                  onChange={handleStudyChange}
                />
              </Col>
              <Col lg="2">
                <Select
                  id="csPubCancerType"
                  label="Cancer Type"
                  value={pubCancerType}
                  options={pubCancerTypeOptions}
                  onChange={handleCancerChange}
                />
              </Col>
              <Col lg="2" />
              <Col lg="2" className="d-flex justify-content-end">
                <Button
                  className="mt-auto mb-3"
                  variant="primary"
                  onClick={() =>
                    calculateR('pub', 'cosineSimilarityPublic', {
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
          </Form>

          <div id="pubPlot">
            {pubErr && (
              <div className="p-3">
                <p>An error has occured. Please verify your input.</p>
              </div>
            )}
            {pubPlotURL && (
              <>
                <hr />
                <Plot
                  className="p-3"
                  plotName={pubPlotURL.split('/').slice(-1)[0]}
                  plotURL={pubPlotURL}
                  txtPath={projectID + pubTxtPath}
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
        onSelect={(tab) => dispatchCosineSimilarity({ display: tab })}
      >
        <Nav variant="tabs">
          <Item>
            <Link eventKey="within" as="button" className="outline-none">
              <strong>CS Within Samples</strong>
            </Link>
          </Item>
          <Item>
            <Link eventKey="reference" as="button" className="outline-none">
              <strong>CS to Reference Signatures</strong>
            </Link>
          </Item>
          {source == 'user' && (
            <Item>
              <Link eventKey="public" as="button" className="outline-none">
                <strong>CS to Public Data</strong>
              </Link>
            </Item>
          )}
        </Nav>
        <Content
          className={`bg-white tab-pane-bordered rounded-0 d-block`}
          style={{ overflowX: 'auto' }}
        >
          {tabs.map(({ key, component }) => (
            <Pane key={key} eventKey={key} className="border-0">
              {component}
            </Pane>
          ))}
        </Content>
      </Container>
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
