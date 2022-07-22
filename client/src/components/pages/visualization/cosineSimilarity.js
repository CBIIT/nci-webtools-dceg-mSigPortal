import React, { useState, useEffect } from 'react';
import { Form, Row, Col, Button, Tab, Nav } from 'react-bootstrap';
import { NavHashLink } from 'react-router-hash-link';
// import { useForm } from 'react-hook-form';
import { useVisualizationOptionsQuery } from '../../../services/store/rootApi';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import SvgContainer from '../../controls/svgContainer/svgContainer';
import CustomSelect from '../../controls/select/select-old';
import Description from '../../controls/description/description';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../services/store/visualization';
import { actions as modalActions } from '../../../services/store/modal';
import { defaultMatrix } from '../../../services/utils';

const actions = { ...visualizationActions, ...modalActions };

const { Container, Content, Pane } = Tab;
const { Item, Link } = Nav;

export default function CosineSimilarity({ submitR, getRefSigOptions }) {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.visualization);

  const mergeState = (state) =>
    dispatch(actions.mergeVisualization({ cosineSimilarity: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const { data, error, isFetching } = useVisualizationOptionsQuery();

  const studyOptions = data
    ? // ? Object.keys(data).map((e) => ({ label: e, value: e }))
      Object.keys(data)
    : [];

  const { study, cancer, strategy } = store.publicForm;

  const { source, matrixList, svgList, profileOptions } = store.main;

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
  } = store.cosineSimilarity;

  // const { control, watch } = useForm();

  // check for multiple sample input and disable params if true
  const [multiSample, setMultiSample] = useState(false);
  useEffect(() => {
    if (Object.keys(svgList).length) {
      if (source == 'user') {
        const samples = [
          ...new Set(
            svgList.map((row) => {
              if (row.Filter != 'NA') return `${row.sample}@${row.Filter}`;
              else return row.sample;
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
    mergeState({ [`${type}SubmitOverlay`]: display });
  }

  // get Signature Reference Sets for dropdown options
  async function getSignatureSet(profileType) {
    if (profileType) {
      mergeState({ refSubmitOverlay: true });
      try {
        const response = await getRefSigOptions(profileType);

        if (response.ok) {
          const { output: refSignatureSetOptions } = await response.json();

          mergeState({
            refSignatureSetOptions: refSignatureSetOptions,
            refSignatureSet: refSignatureSetOptions[0],
            refSubmitOverlay: false,
          });
        } else {
          mergeError(await response.json());
          mergeState({ refSubmitOverlay: false });
        }
      } catch (err) {
        mergeError(err.message);
        mergeState({ refSubmitOverlay: false });
      }
    }
  }

  function setOverlay(type, status) {
    mergeState({ [`${type}SubmitOverlay`]: status });
  }

  async function calculateR(type, fn, args) {
    try {
      setOverlay(type, true);

      mergeState({
        [`${type}Err`]: false,
        [`${type}PlotPath`]: '',
        [`${type}TxtPath`]: '',
      });

      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();
        mergeState({ [`${type}Err`]: err });
      } else {
        const { output } = await response.json();
        const { plotPath, txtPath, error, uncaughtError } = output;

        if (plotPath) {
          mergeState({
            [`${type}PlotPath`]: plotPath,
            [`${type}TxtPath`]: txtPath,
          });
        } else {
          mergeState({
            [`${type}Err`]: error || uncaughtError || true,
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
    const withinMatrixOptions = [
      ...new Set(
        svgList
          .filter((e) => e.profileType == profileType)
          .map((e) => e.matrixSize)
      ),
    ].sort((a, b) => a - b);

    mergeState({
      withinProfileType: profileType,
      withinMatrixSize: defaultMatrix(profileType, withinMatrixOptions),
      withinMatrixOptions: withinMatrixOptions,
    });
  }

  function handlePublicProfileType(profileType) {
    const userMatrixOptions = [
      ...new Set(
        matrixList
          .filter((e) => e.profileType == profileType)
          .map((e) => e.matrixSize)
      ),
    ].sort((a, b) => a - b);

    mergeState({
      userProfileType: profileType,
      userMatrixSize: defaultMatrix(profileType, userMatrixOptions),
      userMatrixOptions: userMatrixOptions,
    });
  }

  function handleStudyChange(study) {
    const cancerTypeOptions = Object.keys(data[study]);

    mergeState({
      pubStudy: study,
      pubCancerType: cancerTypeOptions[0],
      pubCancerTypeOptions: cancerTypeOptions,
    });
  }

  function handleCancerChange(cancer) {
    mergeState({
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
              <Col lg="auto">
                <CustomSelect
                  disabled={!multiSample}
                  id="csProfileType"
                  label="Profile Type"
                  value={withinProfileType}
                  options={profileOptions}
                  onChange={handleWithinProfileType}
                />
              </Col>
              <Col lg="auto">
                <CustomSelect
                  disabled={!multiSample}
                  id="csMatrixSize"
                  label="Matrix Size"
                  value={withinMatrixSize}
                  options={withinMatrixOptions}
                  onChange={(matrix) =>
                    mergeState({
                      withinMatrixSize: matrix,
                    })
                  }
                />
              </Col>
              <Col lg="auto" className="d-flex align-bottom">
                <Button
                  className="mt-auto mb-3"
                  disabled={!multiSample}
                  variant="primary"
                  onClick={() => {
                    if (source == 'user') {
                      calculateR('within', 'cosineSimilarityWithin', {
                        matrixFile: matrixList.filter(
                          (e) =>
                            e.profileType == withinProfileType &&
                            e.matrixSize == withinMatrixSize
                        )[0].Path,
                      });
                    } else {
                      calculateR('within', 'cosineSimilarityWithinPublic', {
                        profileType: withinProfileType,
                        matrixSize: withinMatrixSize,
                        study: study.value,
                        cancerType: cancer.value,
                        experimentalStrategy: strategy.value,
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
              <>
                <hr />
                <p className="p-3">
                  An error has occured. Please verify your input.
                </p>
              </>
            )}
            {withinPlotPath && (
              <>
                <hr />
                <SvgContainer
                  className="p-3"
                  downloadName={withinPlotPath.split('/').slice(-1)[0]}
                  plotPath={'web/results/' + withinPlotPath}
                  txtPath={`web/results/${withinTxtPath}`}
                />
                <p className="p-3">
                  The heatmap shows pairwise cosine similarity between samples
                  from the selected profile type. On the x-axis and y-axis are
                  the sample names. This analysis will help to highlight the
                  samples within the same cluster that may have similar
                  mutational signatures.
                </p>
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
              <Col lg="auto">
                <CustomSelect
                  id="csRefProfileType"
                  label="Profile Type"
                  value={refProfileType}
                  options={profileOptions}
                  onChange={(refProfileType) => {
                    mergeState({
                      refProfileType: refProfileType,
                    });
                    getSignatureSet(refProfileType);
                  }}
                />
              </Col>
              <Col lg="auto">
                <CustomSelect
                  id="csRefSignatureSet"
                  label="Reference Signature Set"
                  value={refSignatureSet}
                  options={refSignatureSetOptions}
                  onChange={(refSignatureSet) => {
                    mergeState({
                      refSignatureSet: refSignatureSet,
                    });
                  }}
                />
              </Col>
              <Col lg="auto" className="d-flex">
                <Button
                  className="mt-auto mb-3"
                  variant="primary"
                  onClick={() => {
                    if (source == 'user') {
                      calculateR('ref', 'cosineSimilarityRefSig', {
                        profileType: refProfileType,
                        signatureSet: refSignatureSet,
                        matrixFile: matrixList.filter(
                          (row) =>
                            row.profileType == refProfileType &&
                            row.matrixSize ==
                              defaultMatrix(refProfileType, ['96', '78', '83'])
                        )[0].Path,
                      });
                    } else {
                      calculateR('ref', 'cosineSimilarityRefSigPublic', {
                        profileType: refProfileType,
                        signatureSet: refSignatureSet,
                        study: study.value,
                        cancerType: cancer.value,
                        experimentalStrategy: strategy.value,
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
              <>
                <hr />
                <p className="p-3">
                  An error has occured. Please verify your input.
                </p>
              </>
            )}
            {refPlotPath && (
              <>
                <hr />
                <SvgContainer
                  className="p-3"
                  downloadName={refPlotPath.split('/').slice(-1)[0]}
                  plotPath={`web/results/${refPlotPath}`}
                  txtPath={`web/results/${refTxtPath}`}
                />
                <p className="p-3">
                  The following heatmap shows pairwise cosine similarity between
                  the mutational profiles of given samples and the selected
                  reference signature set. On the x-axis and y-axis are the
                  reference signature names and the sample names, respectively.
                  This analysis will identify potential dominant mutational
                  signatures in selected samples.
                </p>
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
              <Col lg="auto">
                <CustomSelect
                  id="csUserProfileType"
                  label="Profile Type"
                  value={userProfileType}
                  options={profileOptions}
                  onChange={handlePublicProfileType}
                />
              </Col>
              <Col lg="auto">
                <CustomSelect
                  id="csUserMatrixSize"
                  label="Matrix Size"
                  value={userMatrixSize}
                  options={userMatrixOptions}
                  onChange={(matrix) =>
                    mergeState({
                      userMatrixSize: matrix,
                    })
                  }
                />
              </Col>
              <Col lg="auto">
                <CustomSelect
                  id="csPubStudy"
                  label="Study"
                  value={pubStudy || study}
                  options={studyOptions}
                  onChange={handleStudyChange}
                />
              </Col>
              <Col lg="auto">
                <CustomSelect
                  id="csPubCancerType"
                  label="Cancer Type or Group"
                  value={pubCancerType || cancer}
                  options={pubCancerTypeOptions}
                  onChange={handleCancerChange}
                />
              </Col>
              <Col lg="auto" className="d-flex">
                <Button
                  className="mt-auto mb-3"
                  variant="primary"
                  onClick={() =>
                    calculateR('pub', 'cosineSimilarityPublic', {
                      matrixFile: matrixList.filter(
                        (row) =>
                          row.profileType == userProfileType &&
                          row.matrixSize == userMatrixSize
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
              <>
                <hr />
                <p className="p-3">
                  An error has occured. Please verify your input.
                </p>
              </>
            )}
            {pubPlotPath && (
              <>
                <hr />
                <SvgContainer
                  className="p-3"
                  downloadName={pubPlotPath.split('/').slice(-1)[0]}
                  plotPath={`web/results/${pubPlotPath}`}
                  txtPath={`web/results/${pubTxtPath}`}
                />
                <p className="p-3">
                  The following heatmap shows pairwise cosine similarity between
                  mutational profiles of samples from the user input and public
                  dataset based on the selected study and cancer type. On the
                  x-axis and y-axis are the sample names from the public dataset
                  and user input, respectively.
                </p>
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
        onSelect={(tab) => mergeState({ display: tab })}
      >
        <Nav variant="tabs">
          <Item>
            <Link eventKey="within" as="button" className="outline-none">
              <strong>CS Between Samples</strong>
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
              <div className="p-3">
                <Description
                  less="Cosine similarity is a measure of the similarity of two
                      matrices, which can be helpful to compare two mutational
                      profiles or signatures."
                  more={
                    <span>
                      Below you can explore cosine similarity between sample
                      profiles (CS Between Samples), cosine similarity between
                      sample profiles and reference signatures (CS to Reference
                      Signatures), or, if using your own data, cosine similarity
                      between profiles from your input data and profiles from
                      public data (CS to Public Data). Simply use the dropdown
                      menus to select a [Profile Type], [Matrix Size], or
                      [Reference Signature Set]. Click here to learn more about
                      cosine similarity. Click{' '}
                      <NavHashLink to="/faq#cosine-similarity">
                        here
                      </NavHashLink>{' '}
                      to learn more about cosine similarity.
                    </span>
                  }
                />
              </div>
              <hr />
              {component}
            </Pane>
          ))}
        </Content>
      </Container>
    </div>
  );
}
