import React, { useState, useEffect } from 'react';
import { Form, Row, Col, Button, Tab, Nav } from 'react-bootstrap';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';
import Select from '../../controls/select/select';
import Description from '../../controls/description/description';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../services/store/visualization';
import { actions as modalActions } from '../../../services/store/modal';
import { defaultMatrix } from '../../../services/utils';

const actions = { ...visualizationActions, ...modalActions };
const { Container, Content, Pane } = Tab;
const { Item, Link } = Nav;

export default function PCA({ submitR, getRefSigOptions }) {
  const dispatch = useDispatch();
  const visualization = useSelector((state) => state.visualization);
  const mergePCA = (state) =>
    dispatch(actions.mergeVisualization({ pca: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const {
    source,
    study,
    studyOptions,
    cancerType,
    pubExperimentalStrategy,
    pDataOptions,
    matrixList,
    projectID,
    svgList,
  } = visualization.state;
  const { profileOptions } = visualization.mutationalProfiles;
  const {
    profileType,
    signatureSet,
    signatureSetOptions,
    pca1,
    pca2,
    pca3,
    heatmap,
    pca2Data,
    pca3Data,
    heatmapData,
    pcaErr,
    debugR,
    submitOverlay,
    userProfileType,
    userMatrixSize,
    userMatrixOptions,
    pubStudy,
    pubCancerType,
    pubCancerTypeOptions,
    pubPca1,
    pubPca2,
    pubPca3,
    pubPca2Data,
    pubPca3Data,
    display,
    pubPcaErr,
    pubSubmitOverlay,
  } = visualization.pca;

  const [multiSample, setMultiSample] = useState(false);

  // check for multiple sample input and disable params if true
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

  function setOverlay(type, status) {
    type == 'within'
      ? mergePCA({ submitOverlay: status })
      : mergePCA({ pubSubmitOverlay: status });
  }

  // get Signature Reference Sets for dropdown options
  async function getSignatureSet(profileType) {
    try {
      if (profileType) {
        setOverlay('within', true);

        const response = await getRefSigOptions(profileType);

        if (response.ok) {
          const { output: signatureSetOptions } = await response.json();

          mergePCA({
            signatureSetOptions: signatureSetOptions,
            signatureSet: signatureSetOptions[0],
          });
        } else {
          mergeError(await response.json());
        }
      }
    } catch (err) {
      mergeError(err.message);
    } finally {
      setOverlay('within', false);
    }
  }

  async function calculateR(type, fn, args) {
    try {
      setOverlay(type, true);
      if (type == 'within') {
        mergePCA({
          debugR: '',
          pcaErr: false,
          pca1: '',
          pca2: '',
          pca3: '',
          heatmap: '',
        });
      } else {
        mergePCA({
          debugR: '',
          pubPcaErr: false,
          pubPca1: '',
          pubPca2: '',
          pubPca3: '',
        });
      }

      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();
        mergePCA({ debugR: err });
      } else {
        const { output } = await response.json();

        if (output.pca1) {
          if (type == 'within') {
            mergePCA({
              pca1: output.pca1,
              pca2: output.pca2,
              pca3: output.pca3,
              heatmap: output.heatmap,
              pca2Data: output.pca2Data,
              pca3Data: output.pca3Data,
              heatmapData: output.heatmapData,
            });
          } else {
            mergePCA({
              pubPca1: output.pca1,
              pubPca2: output.pca2,
              pubPca3: output.pca3,
              pubPca2Data: output.pca2Data,
              pubPca3Data: output.pca3Data,
            });
          }
        } else {
          if (type == 'within') {
            mergePCA({ pcaErr: output.error || output.uncaughtError || true });
          } else {
            mergePCA({
              pubPcaErr: output.error || output.uncaughtError || true,
            });
          }
        }
      }
    } catch (err) {
      mergeError(err.message);
    } finally {
      setOverlay(type, false);
    }
  }

  function handleProfileType(profileType) {
    const userMatrixOptions = [
      ...new Set(
        svgList
          .filter((plot) => plot.Profile_Type == profileType)
          .map((plot) => plot.Matrix_Size)
      ),
    ].sort((a, b) => a - b);

    mergePCA({
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

    mergePCA({
      pubStudy: study,
      pubCancerType: cancerTypeOptions[0],
      pubCancerTypeOptions: cancerTypeOptions,
    });
  }

  function handleCancerChange(cancer) {
    mergePCA({
      pubCancerType: cancer,
    });
  }
  let tabs = [
    {
      key: 'within',
      component: (
        <div>
          <Form className="p-3">
            <LoadingOverlay active={submitOverlay} />
            <Row>
              <Col lg="auto">
                <Select
                  disabled={!multiSample}
                  id="pcaProfileType"
                  label="Profile Type"
                  value={profileType}
                  options={profileOptions}
                  onChange={(profileType) => {
                    mergePCA({
                      profileType: profileType,
                    });
                    getSignatureSet(profileType);
                  }}
                />
              </Col>

              <Col lg="auto">
                <Select
                  disabled={!multiSample}
                  id="pcaRefSet"
                  label="Reference Signature Set"
                  value={signatureSet}
                  options={signatureSetOptions}
                  onChange={(signatureSet) => {
                    mergePCA({
                      signatureSet: signatureSet,
                    });
                  }}
                />
              </Col>
              <Col lg="aitp" className="d-flex">
                <Button
                  className="mt-auto mb-3"
                  disabled={!multiSample}
                  variant="primary"
                  onClick={() => {
                    if (source == 'user') {
                      calculateR('within', 'pca', {
                        profileType: profileType,
                        signatureSet: signatureSet,
                        matrixFile: matrixList.filter(
                          (row) =>
                            row.Profile_Type == profileType &&
                            row.Matrix_Size ==
                              defaultMatrix(profileType, ['96', '78', '83'])
                        )[0].Path,
                      });
                    } else {
                      calculateR('within', 'pcaPublic', {
                        profileType: profileType,
                        signatureSet: signatureSet,
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
          {pcaErr && (
            <p className="p-3">
              An error has occured. Please verify your input.
            </p>
          )}

          {pca1 && (
            <div id="pca1Plot">
              <Plot
                className="p-3"
                downloadName={pca1.split('/').slice(-1)[0]}
                plotPath={'api/results/' + pca1}
              />
              <p className="p-3">
                The bar plot illustrates each of the principal components on the
                x-axis and the percentage of variation that each component
                explains on the y-axis.
              </p>
            </div>
          )}

          {pca2 && (
            <div id="pca2Plot">
              <hr />
              <Plot
                className="p-3"
                downloadName={pca2.split('/').slice(-1)[0]}
                plotPath={'api/results/' + pca2}
                txtPath={`api/results/${pca2Data}`}
              />
              <p className="p-3">
                The individual PCA plot based on the top two principal
                components helps to explain a majority of the variation in
                selected or input data. Each dot on the plot is a sample. The
                legend on the right denotes the percent contribution (contrib)
                of each sample to the principal components on the graph (Dim 1
                and Dim 2).
              </p>
            </div>
          )}

          {pca3 && (
            <div id="pca3Plot">
              <hr />
              <Plot
                className="p-3"
                downloadName={pca3.split('/').slice(-1)[0]}
                plotPath={'api/results/' + pca3}
                txtPath={`api/results/${pca3Data}`}
              />
              <p className="p-3">
                The variable PCA plot based on the top two principal components
                helps to explain a majority of the variation in the data. The
                legend on the right denotes the percent contribution (contrib)
                of each mutation type to the principal components on the graph
                (Dim 1 and Dim 2).
              </p>
            </div>
          )}

          {heatmap && (
            <div id="heatmapPlot">
              <hr />
              <Plot
                className="p-3"
                downloadName={heatmap.split('/').slice(-1)[0]}
                plotPath={'api/results/' + heatmap}
                txtPath={`api/results/${heatmapData}`}
              />
              <p className="p-3">
                The heatmap shows cosine similarity between each principal
                component and each mutational signature in the selected
                reference signature set. Brighter colors denote higher levels of
                cosine similarity between the principal component and the
                mutational signature. Red dots in some of the boxes indicate a
                cosine similarity less than 0, denoting a negative correlation
                (e.g., a box with cosine similarity value of 0.8 and marked with
                a red dot indicates a true cosine similarity of -0.8).
              </p>
            </div>
          )}
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
                <Select
                  id="pcaPubProfile"
                  label="Profile Type"
                  value={userProfileType}
                  options={profileOptions}
                  onChange={handleProfileType}
                />
              </Col>
              <Col lg="auto">
                <Select
                  id="pcaPubMatrixSize"
                  label="Matrix Size"
                  value={userMatrixSize}
                  options={userMatrixOptions}
                  onChange={(matrix) => {
                    mergePCA({ userMatrixSize: matrix });
                  }}
                />
              </Col>
              <Col lg="auto">
                <Select
                  id="pcaPubStudy"
                  label="Study"
                  value={pubStudy}
                  options={studyOptions}
                  onChange={(study) => handleStudyChange(study)}
                />
              </Col>
              <Col lg="auto">
                <Select
                  id="pcaPubCancerType"
                  label="Cancer Type or Group"
                  value={pubCancerType}
                  options={pubCancerTypeOptions}
                  onChange={handleCancerChange}
                />
              </Col>
              <Col lg="auto" className="d-flex">
                <Button
                  className="mt-auto mb-3"
                  variant="primary"
                  onClick={() =>
                    calculateR('pub', 'pcaWithPublic', {
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

          {pubPcaErr && (
            <p className="p-3">
              An error has occured. Please verify your input.
            </p>
          )}

          {pubPca1 && (
            <div id="pubPca1Plot">
              <Plot
                className="p-3"
                downloadName={pubPca1.split('/').slice(-1)[0]}
                plotPath={'api/results/' + pubPca1}
              />
              <p className="p-3">
                The bar plot illustrates each of the principal components on the
                x-axis and the percentage of variation that each component
                explains on the y-axis.
              </p>
            </div>
          )}

          {pubPca2 && (
            <div id="pubPca2Plot">
              <hr />
              <Plot
                className="p-3"
                downloadName={pubPca2.split('/').slice(-1)[0]}
                plotPath={'api/results/' + pubPca2}
                txtPath={`api/results/${pubPca2Data}`}
              />
              <p className="p-3">
                The individual PCA plot based on the top two principal
                components helps to explain a majority of the variation in
                selected or input data. Each dot on the plot is a sample. The
                legend on the right denotes the percent contribution (contrib)
                of each sample to the principal components on the graph (Dim 1
                and Dim 2).
              </p>
            </div>
          )}

          {pubPca3 && (
            <div id="pubPca3Plot">
              <hr />
              <Plot
                className="p-3"
                downloadName={pubPca3.split('/').slice(-1)[0]}
                plotPath={'api/results/' + pubPca3}
                txtPath={`api/results/${pubPca3Data}`}
              />
              <p className="p-3">
                The variable PCA plot based on the top two principal components
                helps to explain a majority of the variation in the data. The
                legend on the right denotes the percent contribution (contrib)
                of each mutation type to the principal components on the graph
                (Dim 1 and Dim 2).
              </p>
            </div>
          )}
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
        onSelect={(tab) => mergePCA({ display: tab })}
      >
        <Nav variant="tabs">
          <Item>
            <Link eventKey="within" as="button" className="outline-none">
              <strong>PCA Between Samples</strong>
            </Link>
          </Item>
          {source == 'user' && (
            <Item>
              <Link eventKey="public" as="button" className="outline-none">
                <strong>PCA with Public Data</strong>
              </Link>
            </Item>
          )}
        </Nav>
        <Content
          className={`bg-white tab-pane-bordered rounded-0 d-block`}
          style={{ overflowX: 'auto', minHeight: '500px' }}
        >
          <div className="p-3">
            <Description
              less="Below you can conduct PCA analysis between samples, or a PCA with Public Data (for user input data only)."
              more="PCA stands for Principal Component Analysis, which helps to explain the variation found in the data through the establishment of different principal components. Each principal component can also be used to compare with known mutational signatures."
            />
          </div>
          <hr />
          {tabs.map(({ key, component }) => (
            <Pane key={key} eventKey={key} className="border-0">
              {component}
            </Pane>
          ))}
        </Content>
      </Container>
    </div>
  );
}
