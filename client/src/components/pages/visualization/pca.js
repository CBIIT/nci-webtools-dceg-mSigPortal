import React, { useEffect } from 'react';
import { Form, Row, Col, Button, Accordion, Card } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import { dispatchError, dispatchPCA } from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';
import Debug from '../../controls/debug/debug';
import Select from '../../controls/select/select';

const { Header, Body } = Card;
const { Toggle, Collapse } = Accordion;

export default function PCA({ submitR, getRefSigOptions }) {
  const { matrixList } = useSelector((state) => state.visualizeResults);
  const {
    source,
    study,
    studyOptions,
    cancerType,
    pubExperimentalStrategy,
    pDataOptions,
  } = useSelector((state) => state.visualize);
  const { profileOptions } = useSelector((state) => state.mutationalProfiles);
  const rootURL = window.location.pathname;
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
    displayPCA,
    pcaErr,
    debugR,
    displayDebug,
    submitOverlay,
    pubPCA,

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
    displayPub,
    pubPcaErr,
    pubSubmitOverlay,
  } = useSelector((state) => state.pca);

  async function setRPlot(plotPath, type) {
    if (plotPath) {
      !type.includes('pub')
        ? dispatchPCA({ submitOverlay: true })
        : dispatchPCA({ pubSubmitOverlay: true });
      try {
        const response = await fetch(`${rootURL}results/${plotPath}`);

        if (!response.ok) {
          // console.log(await response.json());
        } else {
          const pic = await response.blob();
          const objectURL = URL.createObjectURL(pic);

          if (type == 'pca1') {
            if (pca1URL) URL.revokeObjectURL(pca1URL);
            dispatchPCA({ pca1URL: objectURL });
          } else if (type == 'pca2') {
            if (pca2URL) URL.revokeObjectURL(pca2URL);
            dispatchPCA({ pca2URL: objectURL });
          } else if (type == 'pca3') {
            if (pca3URL) URL.revokeObjectURL(pca3URL);
            dispatchPCA({ pca3URL: objectURL });
          } else if (type == 'heatmap') {
            if (heatmapURL) URL.revokeObjectURL(heatmapURL);
            dispatchPCA({ heatmapURL: objectURL });
          } else if (type == 'pubPca1') {
            if (pubPca1URL) URL.revokeObjectURL(pubPca1);
            dispatchPCA({ pubPca1URL: objectURL });
          } else if (type == 'pubPca2') {
            if (pubPca2URL) URL.revokeObjectURL(pubPca2URL);
            dispatchPCA({ pubPca2URL: objectURL });
          } else if (type == 'pubPca3') {
            if (pubPca3URL) URL.revokeObjectURL(pubPca3URL);
            dispatchPCA({ pubPca3URL: objectURL });
          }
        }
      } catch (err) {
        dispatchError(err);
      }
      !type.includes('pub')
        ? dispatchPCA({ submitOverlay: false })
        : dispatchPCA({ pubSubmitOverlay: false });
    }
  }

  // get Signature Reference Sets for dropdown options
  async function getSignatureSet(profileType) {
    if (profileType) {
      dispatchPCA({ submitOverlay: true });
      try {
        const response = await getRefSigOptions(profileType);

        if (response.ok) {
          const signatureSetOptions = await response.json();

          dispatchPCA({
            signatureSetOptions: signatureSetOptions,
            signatureSet: signatureSetOptions[0],
            submitOverlay: false,
          });
        } else {
          dispatchError(await response.json());
          dispatchPCA({ submitOverlay: false });
        }
      } catch (err) {
        dispatchError(err);
        dispatchPCA({ submitOverlay: false });
      }
    }
  }

  async function calculateR(fn, args) {
    if (!fn.includes('With')) {
      dispatchPCA({
        submitOverlay: true,
        debugR: '',
        pcaErr: false,
      });
    } else {
      dispatchPCA({
        pubSubmitOverlay: true,
        debugR: '',
        pubPcaErr: false,
      });
    }

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();

        dispatchPCA({ debugR: err });
        if (!fn.includes('With')) {
          dispatchPCA({ submitOverlay: false });
        } else {
          dispatchPCA({
            pubSubmitOverlay: false,
          });
        }
      } else {
        const { debugR, output } = await response.json();

        if (!fn.includes('With')) {
          if (output.pca1) {
            dispatchPCA({
              debugR: debugR,
              submitOverlay: false,
              pca1: output.pca1,
              pca2: output.pca2,
              pca3: output.pca3,
              heatmap: output.heatmap,
              pca2Data: output.pca2Data,
              pca3Data: output.pca3Data,
              heatmapData: output.heatmapData,
            });
            setRPlot(output.pca1, 'pca1');
            setRPlot(output.pca2, 'pca2');
            setRPlot(output.pca3, 'pca3');
            setRPlot(output.heatmap, 'heatmap');
          } else {
            if (pca1URL) URL.revokeObjectURL(pca1URL);
            if (pca2URL) URL.revokeObjectURL(pca2URL);
            if (pca3URL) URL.revokeObjectURL(pca3URL);
            if (heatmapURL) URL.revokeObjectURL(heatmapURL);
            dispatchPCA({
              debugR: debugR,
              submitOverlay: false,
              pcaErr: true,
              pca1URL: '',
              pca2URL: '',
              pca3URL: '',
              heatmapURL: '',
            });
          }
        } else {
          if (output.pubPca1) {
            dispatchPCA({
              debugR: debugR,
              pubSubmitOverlay: false,
              pubPca1: output.pca1,
              pubPca2: output.pca2,
              pubPca3: output.pca3,
              pubPca2Data: output.pca2Data,
              pubPca3Data: output.pca3Data,
            });
            setRPlot(output.pca1, 'pubPca1');
            setRPlot(output.pca2, 'pubPca2');
            setRPlot(output.pca3, 'pubPca3');
          } else {
            if (pubPca1URL) URL.revokeObjectURL(pubPca1URL);
            if (pubPca2URL) URL.revokeObjectURL(pubPca2URL);
            if (pubPca3URL) URL.revokeObjectURL(pubPca3URL);
            dispatchPCA({
              debugR: debugR,
              pubPcaErr: true,
              pubSubmitOverlay: false,
              pubPca1URL: '',
              pubPca2URL: '',
              pubPca3URL: '',
            });
          }
        }
      }
    } catch (err) {
      dispatchError(err);
      if (!fn.includes('With')) {
        dispatchPCA({ submitOverlay: false });
      } else {
        dispatchPCA({
          pubSubmitOverlay: false,
        });
      }
    }
  }

  function handleProfileType(profileType) {
    const userMatrixOptions = [
      ...new Set(
        matrixList
          .filter((matrix) => matrix.Profile_Type == profileType)
          .map((matrix) => matrix.Matrix_Size)
      ),
    ];

    dispatchPCA({
      userProfileType: profileType,
      userMatrixSize: userMatrixOptions[0],
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

    dispatchPCA({
      pubStudy: study,
      pubCancerType: cancerTypeOptions[0],
      pubCancerTypeOptions: cancerTypeOptions,
    });
  }

  function handleCancerChange(cancer) {
    dispatchPCA({
      pubCancerType: cancer,
    });
  }

  return (
    <div>
      <Accordion defaultActiveKey="0">
        <Card>
          <Toggle
            className="font-weight-bold"
            as={Header}
            eventKey="0"
            onClick={() =>
              dispatchPCA({
                displayPCA: !displayPCA,
              })
            }
          >
            {displayPCA == true ? (
              <FontAwesomeIcon icon={faMinus} />
            ) : (
              <FontAwesomeIcon icon={faPlus} />
            )}{' '}
            PCA Within Samples
          </Toggle>
          <Collapse eventKey="0">
            <Body>
              <Form>
                <LoadingOverlay active={submitOverlay} />
                <div>
                  <Row className="justify-content-center">
                    <Col sm="5">
                      <Select
                        id="pcaProfileType"
                        label="Profile Type"
                        value={profileType}
                        options={profileOptions}
                        onChange={(profileType) => {
                          dispatchPCA({
                            profileType: profileType,
                          });
                          getSignatureSet(profileType);
                        }}
                      />
                    </Col>

                    <Col sm="5">
                      <Select
                        id="pcaRefSet"
                        label="Reference Signature Set"
                        value={signatureSet}
                        options={signatureSetOptions}
                        onChange={(signatureSet) => {
                          dispatchPCA({
                            signatureSet: signatureSet,
                          });
                        }}
                      />
                    </Col>
                    <Col sm="1" className="m-auto">
                      <Button
                        variant="primary"
                        onClick={() => {
                          if (source == 'user') {
                            calculateR('pca', {
                              profileType: profileType,
                              signatureSet: signatureSet,
                              matrixList: JSON.stringify(
                                matrixList.filter(
                                  (matrix) => matrix.Profile_Type == profileType
                                )
                              ),
                            });
                          } else {
                            calculateR('pcaPublic', {
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

                  <div id="pca1Plot">
                    <div style={{ display: pcaErr ? 'block' : 'none' }}>
                      <p>
                        An error has occured. Check the debug section for more
                        info.
                      </p>
                    </div>
                    <div
                      className="my-4"
                      style={{ display: pca1URL ? 'block' : 'none' }}
                    >
                      <Plot
                        plotName={pca1.split('/').slice(-1)[0]}
                        plotURL={pca1URL}
                      />
                    </div>
                  </div>

                  <div id="pca2Plot">
                    <div
                      className="my-4"
                      style={{ display: pca2URL ? 'block' : 'none' }}
                    >
                      <Plot
                        plotName={pca2.split('/').slice(-1)[0]}
                        plotURL={pca2URL}
                        txtPath={pca2Data}
                      />
                    </div>
                  </div>

                  <div id="pca3Plot">
                    <div
                      className="my-4"
                      style={{ display: pca3URL ? 'block' : 'none' }}
                    >
                      <Plot
                        plotName={pca3.split('/').slice(-1)[0]}
                        plotURL={pca3URL}
                        txtPath={pca3Data}
                      />
                    </div>
                  </div>

                  <div id="heatmapPlot">
                    <div
                      className="my-4"
                      style={{ display: heatmapURL ? 'block' : 'none' }}
                    >
                      <Plot
                        plotName={heatmap.split('/').slice(-1)[0]}
                        plotURL={heatmapURL}
                        txtPath={heatmapData}
                      />
                    </div>
                  </div>
                </div>
              </Form>
            </Body>
          </Collapse>
        </Card>
      </Accordion>

      {source == 'user' && (
        <Accordion defaultActiveKey="1">
          <Card>
            <Toggle
              className="font-weight-bold"
              as={Header}
              eventKey="1"
              onClick={() => dispatchPCA({ displayPub: !displayPub })}
            >
              {displayPub == true ? (
                <FontAwesomeIcon icon={faMinus} />
              ) : (
                <FontAwesomeIcon icon={faPlus} />
              )}{' '}
              PCA with Public Data
            </Toggle>
            <Collapse eventKey="1">
              <Body>
                <Form>
                  <LoadingOverlay active={pubSubmitOverlay} />
                  <div>
                    <Row className="justify-content-center">
                      <Col sm="2">
                        <Select
                          id="pcaPubProfile"
                          label="Profile Type"
                          value={userProfileType}
                          options={profileOptions}
                          onChange={handleProfileType}
                        />
                      </Col>
                      <Col sm="2">
                        <Select
                          id="pcaPubMatrixSize"
                          label="Matrix Size"
                          value={userMatrixSize}
                          options={userMatrixOptions}
                          onChange={(matrix) => {
                            dispatchPCA({ userMatrixSize: matrix });
                          }}
                        />
                      </Col>
                      <Col sm="2">
                        <Select
                          id="pcaPubStudy"
                          label="Study"
                          value={pubStudy}
                          options={studyOptions}
                          onChange={(study) => handleStudyChange(study)}
                        />
                      </Col>
                      <Col sm="4">
                        <Select
                          id="pcaPubCancerType"
                          label="Cancer Type"
                          value={pubCancerType}
                          options={pubCancerTypeOptions}
                          onChange={handleCancerChange}
                        />
                      </Col>
                      <Col sm="1" className="m-auto">
                        <Button
                          variant="primary"
                          onClick={() =>
                            calculateR('pcaWithPublic', {
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

                    <div id="pubPca1Plot">
                      <div style={{ display: pubPcaErr ? 'block' : 'none' }}>
                        <p>
                          An error has occured. Check the debug section for more
                          info.
                        </p>
                      </div>
                      <div
                        className="my-4"
                        style={{ display: pubPca1URL ? 'block' : 'none' }}
                      >
                        <Plot
                          plotName={pubPca1.split('/').slice(-1)[0]}
                          plotURL={pubPca1URL}
                        />
                      </div>
                    </div>

                    <div id="pubPca2Plot">
                      <div
                        className="my-4"
                        style={{ display: pubPca2URL ? 'block' : 'none' }}
                      >
                        <Plot
                          plotName={pubPca2.split('/').slice(-1)[0]}
                          plotURL={pubPca2URL}
                          txtPath={pubPca2Data}
                        />
                      </div>
                    </div>

                    <div id="pubPca3Plot">
                      <div
                        className="my-4"
                        style={{ display: pubPca3URL ? 'block' : 'none' }}
                      >
                        <Plot
                          plotName={pubPca3.split('/').slice(-1)[0]}
                          plotURL={pubPca3URL}
                          txtPath={pubPca3Data}
                        />
                      </div>
                    </div>
                  </div>
                </Form>
              </Body>
            </Collapse>
          </Card>
        </Accordion>
      )}

      <Debug msg={debugR} />
    </div>
  );
}
