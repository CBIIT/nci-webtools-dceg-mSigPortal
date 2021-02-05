import React, { useState, useEffect } from 'react';
import { Form, Row, Col, Button, Tab, Nav } from 'react-bootstrap';
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
const { Container, Content, Pane } = Tab;
const { Item, Link } = Nav;

export default function PCA({ submitR, getRefSigOptions }) {
  const dispatch = useDispatch();
  const visualization = useSelector((state) => state.visualization);
  const mergePCA = (state) =>
    dispatch(actions.mergeVisualization({ pca: state }));
  const mergeError = (state) => dispatch(actions.mergeModal({ error: state }));

  const { matrixList, projectID, svgList } = visualization.results;
  const {
    source,
    study,
    studyOptions,
    cancerType,
    pubExperimentalStrategy,
    pDataOptions,
  } = visualization.visualize;
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
    pca1URL,
    pca2URL,
    pca3URL,
    heatmapURL,
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
    pubPca1URL,
    pubPca2URL,
    pubPca3URL,
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
            svgList.data.map((row) => {
              if (value2d(row, 'Filter', svgList.columns) != 'NA')
                return `${value2d(
                  row,
                  'Sample_Name',
                  svgList.columns
                )}@${value2d(row, 'Filter', svgList.columns)}`;
              else return value2d(row, 'Sample_Name', svgList.columns);
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
    pca1 ? setRPlot(pca1, 'pca1') : clearPlot('pca1');
    pca2 ? setRPlot(pca2, 'pca2') : clearPlot('pca2');
    pca3 ? setRPlot(pca3, 'pca3') : clearPlot('pca3');
    heatmap ? setRPlot(heatmap, 'heatmap') : clearPlot('heatmap');
    pubPca1 ? setRPlot(pubPca1, 'pubPca1') : clearPlot('pubPca1');
    pubPca2 ? setRPlot(pubPca2, 'pubPca2') : clearPlot('pubPca2');
    pubPca3 ? setRPlot(pubPca3, 'pubPca3') : clearPlot('pubPca3');
  }, [pca1, pca2, pca3, heatmap, pubPca1, pubPca2, pubPca3]);

  function setOverlay(type, status) {
    type == 'within'
      ? mergePCA({ submitOverlay: status })
      : mergePCA({ pubSubmitOverlay: status });
  }

  async function setRPlot(plotPath, type) {
    try {
      const response = await fetch(`api/results/${projectID}${plotPath}`);

      if (!response.ok) {
        // console.log(await response.json());
      } else {
        const pic = await response.blob();
        const objectURL = URL.createObjectURL(pic);

        if (visualization.pca[`${type}URL`]) URL.revokeObjectURL(visualization.pca[`${type}URL`]);
        mergePCA({ [`${type}URL`]: objectURL });
      }
    } catch (err) {
      mergeError({ visible: true, message: err.message });
;
    }
  }

  function clearPlot(type) {
    URL.revokeObjectURL(visualization.pca[`${type}URL`]);
    mergePCA({ [`${type}URL`]: '' });
  }

  // get Signature Reference Sets for dropdown options
  async function getSignatureSet(profileType) {
    try {
      if (profileType) {
        setOverlay('within', true);

        const response = await getRefSigOptions(profileType);

        if (response.ok) {
          const signatureSetOptions = await response.json();

          mergePCA({
            signatureSetOptions: signatureSetOptions,
            signatureSet: signatureSetOptions[0],
          });
        } else {
          mergeError(await response.json());
        }
      }
    } catch (err) {
      mergeError({ visible: true, message: err.message });
;
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
        const { debugR, output } = await response.json();

        mergePCA({ debugR: debugR });

        if (Object.keys(output).length) {
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
            mergePCA({
              debugR: debugR,
              pcaErr: true,
            });
          } else {
            mergePCA({
              debugR: debugR,
              pubPcaErr: true,
            });
          }
        }
      }
    } catch (err) {
      mergeError({ visible: true, message: err.message });
;
    } finally {
      setOverlay(type, false);
    }
  }

  function handleProfileType(profileType) {
    const userMatrixOptions = unique2d(
      'Matrix_Size',
      matrixList.columns,
      filter2d(profileType, matrixList.data)
    );

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
              <Col lg="2">
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

              <Col lg="3">
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
              <Col lg="7" className="d-flex justify-content-end">
                <Button
                  className="mt-auto mb-3"
                  disabled={!multiSample}
                  variant="primary"
                  onClick={() => {
                    if (source == 'user') {
                      calculateR('within', 'pca', {
                        profileType: profileType,
                        signatureSet: signatureSet,
                        matrixFile: value2d(
                          filter2d(
                            [
                              profileType,
                              defaultMatrix(profileType, ['96', '78', '83']),
                            ],
                            matrixList.data
                          )[0],
                          'Path',
                          matrixList.columns
                        ),
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

          {pcaErr && (
            <div>
              <hr />
              <p className="p-3">
                An error has occured. Please verify your input.
              </p>
            </div>
          )}

          {pca1URL && (
            <div id="pca1Plot">
              <hr />
              <Plot
                className="p-3"
                plotName={pca1.split('/').slice(-1)[0]}
                plotURL={pca1URL}
              />
            </div>
          )}

          {pca2URL && (
            <div id="pca2Plot">
              <hr />
              <Plot
                className="p-3"
                plotName={pca2.split('/').slice(-1)[0]}
                plotURL={pca2URL}
                txtPath={projectID + pca2Data}
              />
            </div>
          )}

          {pca3URL && (
            <div id="pca3Plot">
              <hr />
              <Plot
                className="p-3"
                plotName={pca3.split('/').slice(-1)[0]}
                plotURL={pca3URL}
                txtPath={projectID + pca3Data}
              />
            </div>
          )}

          {heatmapURL && (
            <div id="heatmapPlot">
              <hr />
              <Plot
                className="p-3"
                plotName={heatmap.split('/').slice(-1)[0]}
                plotURL={heatmapURL}
                txtPath={projectID + heatmapData}
              />
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
              <Col lg="2">
                <Select
                  id="pcaPubProfile"
                  label="Profile Type"
                  value={userProfileType}
                  options={profileOptions}
                  onChange={handleProfileType}
                />
              </Col>
              <Col lg="2">
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
              <Col lg="2">
                <Select
                  id="pcaPubStudy"
                  label="Study"
                  value={pubStudy}
                  options={studyOptions}
                  onChange={(study) => handleStudyChange(study)}
                />
              </Col>
              <Col lg="3">
                <Select
                  id="pcaPubCancerType"
                  label="Cancer Type"
                  value={pubCancerType}
                  options={pubCancerTypeOptions}
                  onChange={handleCancerChange}
                />
              </Col>
              <Col lg="3" className="d-flex justify-content-end">
                <Button
                  className="mt-auto mb-3"
                  variant="primary"
                  onClick={() =>
                    calculateR('pub', 'pcaWithPublic', {
                      matrixFile: value2d(
                        filter2d(
                          [userProfileType, userMatrixSize],
                          matrixList.data
                        )[0],
                        'Path',
                        matrixList.columns
                      ),
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

          {pubPcaErr && (
            <div>
              <hr />
              <p className="p-3">
                An error has occured. Please verify your input.
              </p>
            </div>
          )}

          {pubPca1URL && (
            <div id="pubPca1Plot">
              <hr />
              <Plot
                className="p-3"
                plotName={pubPca1.split('/').slice(-1)[0]}
                plotURL={pubPca1URL}
              />
            </div>
          )}

          {pubPca2URL && (
            <div id="pubPca2Plot">
              <hr />
              <Plot
                className="p-3"
                plotName={pubPca2.split('/').slice(-1)[0]}
                plotURL={pubPca2URL}
                txtPath={projectID + pubPca2Data}
              />
            </div>
          )}

          {pubPca3URL && (
            <div id="pubPca3Plot">
              <hr />
              <Plot
                className="p-3"
                plotName={pubPca3.split('/').slice(-1)[0]}
                plotURL={pubPca3URL}
                txtPath={projectID + pubPca3Data}
              />
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
              <strong>PCA Within Samples</strong>
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
