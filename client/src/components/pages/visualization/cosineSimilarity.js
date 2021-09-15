import React, { useState, useEffect } from 'react';
import { Form, Row, Col, Button, Tab, Nav } from 'react-bootstrap';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';
import Select from '../../controls/select/select';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../services/store/visualization';
import { actions as modalActions } from '../../../services/store/modal';
import { defaultMatrix } from '../../../services/utils';

const actions = { ...visualizationActions, ...modalActions };

const { Container, Content, Pane } = Tab;
const { Item, Link } = Nav;

export default function CosineSimilarity({ submitR, getRefSigOptions }) {
  const dispatch = useDispatch();
  const visualization = useSelector((state) => state.visualization);

  const mergeCosineSimilarity = (state) =>
    dispatch(actions.mergeVisualization({ cosineSimilarity: state }));
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
  } = visualization.state;
  const { profileOptions } = visualization.mutationalProfiles;

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
    display,
    withinErr,
    refErr,
    pubErr,
    withinSubmitOverlay,
    refSubmitOverlay,
    pubSubmitOverlay,
    debugR,
  } = visualization.cosineSimilarity;

  // check for multiple sample input and disable params if true
  const [multiSample, setMultiSample] = useState(false);
  useEffect(() => {
    if (Object.keys(svgList).length) {
      if (source == 'user') {
        const samples = [
          ...new Set(
            svgList.map((row) => {
              if (row.Filter != 'NA') return `${row.Sample_Name}@${row.Filter}`;
              else return row.Sample_Name;
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

  function setOverlay(type, display) {
    mergeCosineSimilarity({ [`${type}SubmitOverlay`]: display });
  }

  // get Signature Reference Sets for dropdown options
  async function getSignatureSet(profileType) {
    if (profileType) {
      mergeCosineSimilarity({ refSubmitOverlay: true });
      try {
        const response = await getRefSigOptions(profileType);

        if (response.ok) {
          const refSignatureSetOptions = await response.json();

          mergeCosineSimilarity({
            refSignatureSetOptions: refSignatureSetOptions,
            refSignatureSet: refSignatureSetOptions[0],
            refSubmitOverlay: false,
          });
        } else {
          mergeError(await response.json());
          mergeCosineSimilarity({ refSubmitOverlay: false });
        }
      } catch (err) {
        mergeError(err.message);
        mergeCosineSimilarity({ refSubmitOverlay: false });
      }
    }
  }

  function setOverlay(type, status) {
    mergeCosineSimilarity({ [`${type}SubmitOverlay`]: status });
  }

  async function calculateR(type, fn, args) {
    try {
      setOverlay(type, true);

      mergeCosineSimilarity({
        [`${type}Err`]: false,
        [`${type}PlotPath`]: '',
        [`${type}TxtPath`]: '',
        debugR: '',
      });

      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();
        mergeCosineSimilarity({ debugR: err });
      } else {
        const { debugR, output, error } = await response.json();

        if (Object.keys(output).length) {
          mergeCosineSimilarity({
            [`${type}PlotPath`]: output.plotPath,
            [`${type}TxtPath`]: output.txtPath,
          });
        } else {
          mergeCosineSimilarity({
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

  function handleWithinProfileType(profileType) {
    if (source == 'user') {
      const withinMatrixOptions = [
        ...new Set(
          matrixList
            .filter((row) => row.Profile_Type == profileType)
            .map((row) => row.Matrix_Size)
        ),
      ].sort((a, b) => a - b);

      mergeCosineSimilarity({
        withinProfileType: profileType,
        withinMatrixSize: defaultMatrix(profileType, withinMatrixOptions),
        withinMatrixOptions: withinMatrixOptions,
      });
    } else {
      const withinMatrixOptions = [
        ...new Set(
          svgList
            .filter((row) => row.Profile.includes(profileType))
            .map((row) => row.Profile.match(/\d+/gi)[0])
        ),
      ].sort((a, b) => a - b);

      mergeCosineSimilarity({
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
          .filter((row) => row.Profile_Type == profileType)
          .map((row) => row.Matrix_Size)
      ),
    ].sort((a, b) => a - b);

    mergeCosineSimilarity({
      userProfileType: profileType,
      userMatrixSize: defaultMatrix(profileType, userMatrixOptions),
      userMatrixOptions: userMatrixOptions,
    });
  }

  function handleStudyChange(study) {
    const cancerTypeOptions = [
      ...new Set(
        pDataOptions
          .filter((row) => row.Study == study)
          .map((row) => row.Cancer_Type)
      ),
    ];

    mergeCosineSimilarity({
      pubStudy: study,
      pubCancerType: cancerTypeOptions[0],
      pubCancerTypeOptions: cancerTypeOptions,
    });
  }

  function handleCancerChange(cancer) {
    mergeCosineSimilarity({
      pubCancerType: cancer,
    });
  }

  let tabs = [
    {
      key: 'within',
      component: (
        <div>
          <div className="p-3">
            <p>
              Cosine similarity is a measure of the similarity of two matrix,
              which can be helpful to compare two mutational profiles or
              signatures. Below you can observe cosine similarity within sample
              profiles (CS Within Samples), cosine similarity between sample
              profiles and reference signatures (CS to Reference Signatures), or
              if using your own data, cosine similarity between profiles from
              your input data and profiles from public data (CS to Public Data).
              Simply use the dropdown menus to select a “Profile Type”, and the
              “Matrix Size”, Reference Signature Set”, or “Study”, depending on
              the calculation being made.
            </p>
            <p>
              Click <a href="#faq">here</a> to learn more about cosine
              similarity.
            </p>
          </div>
          <hr />
          <Form className="p-3">
            <LoadingOverlay active={withinSubmitOverlay} />
            <Row>
              <Col lg="auto">
                <Select
                  disabled={!multiSample}
                  id="csProfileType"
                  label="Profile Type"
                  value={withinProfileType}
                  options={profileOptions}
                  onChange={handleWithinProfileType}
                />
              </Col>
              <Col lg="auto">
                <Select
                  disabled={!multiSample}
                  id="csMatrixSize"
                  label="Matrix Size"
                  value={withinMatrixSize}
                  options={withinMatrixOptions}
                  onChange={(matrix) =>
                    mergeCosineSimilarity({
                      withinMatrixSize: matrix,
                    })
                  }
                />
              </Col>
              <Col lg="auto" className="d-flex align-bottom">
                <Button
                  className="ml-auto mb-auto"
                  disabled={!multiSample}
                  variant="primary"
                  onClick={() => {
                    if (source == 'user') {
                      calculateR('within', 'cosineSimilarityWithin', {
                        matrixFile: matrixList.filter(
                          (row) =>
                            row.Profile_Type == withinProfileType &&
                            row.Matrix_Size == withinMatrixSize
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
          <hr />
          <div id="withinPlot">
            {withinErr && (
              <p className="p-3">
                An error has occured. Please verify your input.
              </p>
            )}
            {withinPlotPath && (
              <Plot
                className="p-3"
                downloadName={withinPlotPath.split('/').slice(-1)[0]}
                plotPath={'api/results/' + projectID + withinPlotPath}
                txtPath={projectID + withinTxtPath}
              />
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
              The CS to Reference Signatures plot Signature plot highlights
              cosine similarity between profile in a given sample and a
              reference signature set. Along the bottom of the heatmap are the
              signatures within the reference signature set selected. Along the
              side of the heatmap are the sample profiled from the input data
              file.
            </p>
          </div>
          <hr />
          <Form className="p-3">
            <LoadingOverlay active={refSubmitOverlay} />
            <Row>
              <Col lg="auto">
                <Select
                  id="csRefProfileType"
                  label="Profile Type"
                  value={refProfileType}
                  options={profileOptions}
                  onChange={(refProfileType) => {
                    mergeCosineSimilarity({
                      refProfileType: refProfileType,
                    });
                    getSignatureSet(refProfileType);
                  }}
                />
              </Col>
              <Col lg="auto">
                <Select
                  id="csRefSignatureSet"
                  label="Reference Signature Set"
                  value={refSignatureSet}
                  options={refSignatureSetOptions}
                  onChange={(refSignatureSet) => {
                    mergeCosineSimilarity({
                      refSignatureSet: refSignatureSet,
                    });
                  }}
                />
              </Col>
              <Col lg="auto" className="d-flex">
                <Button
                  className="ml-auto mb-auto"
                  variant="primary"
                  onClick={() => {
                    if (source == 'user') {
                      calculateR('ref', 'cosineSimilarityRefSig', {
                        profileType: refProfileType,
                        signatureSet: refSignatureSet,
                        matrixFile: matrixList.filter(
                          (row) =>
                            row.Profile_Type == refProfileType &&
                            row.Matrix_Size ==
                              defaultMatrix(refProfileType, ['96', '78', '83'])
                        )[0].Path,
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
          <hr />
          <div id="refPlot">
            {refErr && (
              <p className="p-3">
                An error has occured. Please verify your input.
              </p>
            )}
            {refPlotPath && (
              <Plot
                className="p-3"
                downloadName={refPlotPath.split('/').slice(-1)[0]}
                plotPath={`api/results/${projectID}${refPlotPath}`}
                txtPath={projectID + refTxtPath}
              />
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
              The CS to Public Data plot highlights cosine similarity between
              sample profiles from the user dataset and the public data
              available from the given study and cancer type. Along the bottom
              of the heatmap are the samples from the selected public data.
              Along the side of the heatmap are the samples from the input data
              file.
            </p>
          </div>
          <hr />
          <Form className="p-3">
            <LoadingOverlay active={pubSubmitOverlay} />
            <Row>
              <Col lg="auto">
                <Select
                  id="csUserProfileType"
                  label="Profile Type"
                  value={userProfileType}
                  options={profileOptions}
                  onChange={handlePublicProfileType}
                />
              </Col>
              <Col lg="auto">
                <Select
                  id="csUserMatrixSize"
                  label="Matrix Size"
                  value={userMatrixSize}
                  options={userMatrixOptions}
                  onChange={(matrix) =>
                    mergeCosineSimilarity({
                      userMatrixSize: matrix,
                    })
                  }
                />
              </Col>
              <Col lg="auto">
                <Select
                  id="csPubStudy"
                  label="Study"
                  value={pubStudy || visualization.study}
                  options={studyOptions}
                  onChange={handleStudyChange}
                />
              </Col>
              <Col lg="auto">
                <Select
                  id="csPubCancerType"
                  label="Cancer Type"
                  value={pubCancerType || visualization.cancerType}
                  options={pubCancerTypeOptions}
                  onChange={handleCancerChange}
                />
              </Col>
              <Col lg="auto" className="d-flex">
                <Button
                  className="ml-auto mb-auto"
                  variant="primary"
                  onClick={() =>
                    calculateR('pub', 'cosineSimilarityPublic', {
                      matrixFile: matrixList.filter(
                        (row) =>
                          row.Profile_Type == userProfileType &&
                          row.Matrix_Size == userMatrixSize
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
          <hr />
          <div id="pubPlot">
            {pubErr && (
              <p className="p-3">
                An error has occured. Please verify your input.
              </p>
            )}
            {pubPlotPath && (
              <Plot
                className="p-3"
                downloadName={pubPlotPath.split('/').slice(-1)[0]}
                plotPath={`api/results/${projectID}${pubPlotPath}`}
                txtPath={projectID + pubTxtPath}
              />
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
        onSelect={(tab) => mergeCosineSimilarity({ display: tab })}
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
          style={{ overflowX: 'auto', minHeight: '500px' }}
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
