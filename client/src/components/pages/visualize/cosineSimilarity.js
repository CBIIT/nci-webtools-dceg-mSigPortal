import React, { useEffect } from 'react';
import { Form, Row, Col, Button, Accordion, Card } from 'react-bootstrap';
import Select from 'react-select';
import { useSelector } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import {
  dispatchError,
  dispatchCosineSimilarity,
} from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

const { Group, Label, Control } = Form;
const { Header, Body } = Card;
const { Toggle, Collapse } = Accordion;

export default function CosineSimilarity({
  downloadResults,
  submitR,
  getRefSigOptions,
}) {
  const { source, study, cancerType, pExperimentalStrategy } = useSelector(
    (state) => state.visualize
  );
  const { profileOptions } = useSelector((state) => state.mutationalProfiles);
  const { displayTab, matrixList, svgList } = useSelector(
    (state) => state.visualizeResults
  );
  const rootURL = window.location.pathname;
  const {
    withinProfileType,
    withinMatrixSize,
    withinMatrixOptions,
    refProfileType,
    refSignatureSet,
    refSignatureSetOptions,
    withinPlotPath,
    withinTxtPath,
    refPlotPath,
    refTxtPath,
    withinPlotURL,
    refPlotURL,
    displayWithin,
    displayRefSig,
    withinErr,
    refErr,
    withinSubmitOverlay,
    refSubmitOverlay,
    debugR,
    displayDebug,
  } = useSelector((state) => state.cosineSimilarity);

  const selectFix = {
    styles: {
      menuPortal: (base) => ({ ...base, zIndex: 9999 }),
    },
    menuPortalTarget: document.body,
  };

  function setOverlay(type, display) {
    if (type == 'within') {
      dispatchCosineSimilarity({ withinSubmitOverlay: display });
    } else {
      dispatchCosineSimilarity({ refSubmitOverlay: display });
    }
  }

  async function setRPlot(plotPath, type) {
    setOverlay(type, true);
    if (plotPath) {
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
        } else {
          const pic = await response.blob();
          const objectURL = URL.createObjectURL(pic);

          if (type == 'within') {
            if (withinPlotURL) URL.revokeObjectURL(withinPlotURL);
            dispatchCosineSimilarity({
              withinPlotURL: objectURL,
            });
          } else {
            if (refPlotURL) URL.revokeObjectURL(refPlotURL);
            dispatchCosineSimilarity({
              refPlotURL: objectURL,
            });
          }
          setOverlay(type, false);
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (type == 'within') {
        if (withinPlotURL) URL.revokeObjectURL(withinPlotURL);
        dispatchCosineSimilarity({ withinErr: true, withinPlotURL: '' });
      } else {
        if (refPlotURL) URL.revokeObjectURL(refPlotURL);
        dispatchCosineSimilarity({ refErr: true, refPlotURL: '' });
      }
    }
    setOverlay(type, false);
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

  async function calculateR(fn, args) {
    if (
      fn == 'cosineSimilarityWithin' ||
      fn == 'cosineSimilarityWithinPublic'
    ) {
      console.log(fn);
      dispatchCosineSimilarity({
        withinSubmitOverlay: true,
        withinErr: '',
        debugR: '',
      });
    } else {
      dispatchCosineSimilarity({
        refSubmitOverlay: true,
        refErr: '',
        debugR: '',
      });
    }
    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();
        if (
          fn == 'cosineSimilarityWithin' ||
          fn == 'cosineSimilarityWithinPublic'
        ) {
          dispatchCosineSimilarity({
            withinSubmitOverlay: false,
            debugR: err,
          });
        } else {
          dispatchCosineSimilarity({
            refSubmitOverlay: false,
            debugR: err,
          });
        }
      } else {
        const { debugR, output } = await response.json();

        if (
          fn == 'cosineSimilarityWithin' ||
          fn == 'cosineSimilarityWithinPublic'
        ) {
          dispatchCosineSimilarity({
            debugR: debugR,
            withinSubmitOverlay: false,
            withinPlotPath: output.plotPath,
            withinTxtPath: output.txtPath,
          });
          setRPlot(output.plotPath, 'within');
        } else {
          dispatchCosineSimilarity({
            debugR: debugR,
            refSubmitOverlay: false,
            refPlotPath: output.plotPath,
            refTxtPath: output.txtPath,
          });
          setRPlot(output.plotPath, 'refsig');
        }
      }
    } catch (err) {
      dispatchError(err);
      if (
        fn == 'cosineSimilarityWithin' ||
        fn == 'cosineSimilarityWithinPublic'
      ) {
        dispatchCosineSimilarity({ withinSubmitOverlay: false });
      } else {
        dispatchCosineSimilarity({ refSubmitOverlay: false });
      }
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
      ];

      dispatchCosineSimilarity({
        withinProfileType: profileType,
        withinMatrixSize: withinMatrixOptions[0],
        withinMatrixOptions: withinMatrixOptions,
      });
    } else {
      const withinMatrixOptions = [
        ...new Set(
          svgList
            .filter((plot) => plot.Profile.indexOf(profileType) > -1)
            .map((plot) => plot.Profile.match(/\d+/gi)[0])
        ),
      ];

      dispatchCosineSimilarity({
        withinProfileType: profileType,
        withinMatrixSize: withinMatrixOptions[0],
        withinMatrixOptions: withinMatrixOptions,
      });
    }
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
              dispatchCosineSimilarity({
                displayWithin: !displayWithin,
              })
            }
          >
            {displayWithin == true ? (
              <FontAwesomeIcon icon={faMinus} />
            ) : (
              <FontAwesomeIcon icon={faPlus} />
            )}{' '}
            Cosine Similarity Within Samples
          </Toggle>
          <Collapse eventKey="0">
            <Body>
              <Form>
                <LoadingOverlay active={withinSubmitOverlay} />
                <div>
                  <Row className="justify-content-center">
                    <Col sm="5">
                      <Group controlId="withinProfileType">
                        <Label>Profile Type</Label>
                        <Select
                          options={profileOptions}
                          value={[withinProfileType]}
                          onChange={(profile) =>
                            handleWithinProfileType(profile)
                          }
                          getOptionLabel={(option) => option}
                          getOptionValue={(option) => option}
                          {...selectFix}
                        />
                      </Group>
                    </Col>
                    <Col sm="5">
                      <Label>Matrix Size</Label>
                      <Select
                        options={withinMatrixOptions}
                        value={[withinMatrixSize]}
                        onChange={(matrix) =>
                          dispatchCosineSimilarity({
                            withinMatrixSize: matrix,
                          })
                        }
                        getOptionLabel={(option) => option}
                        getOptionValue={(option) => option}
                        {...selectFix}
                      />
                    </Col>
                    <Col sm="2" className="m-auto">
                      <Button
                        variant="primary"
                        onClick={() => {
                          if (source == 'user') {
                            calculateR('cosineSimilarityWithin', {
                              matrixFile: matrixList.filter(
                                (path) =>
                                  path.Profile_Type == withinProfileType &&
                                  path.Matrix_Size == withinMatrixSize
                              )[0].Path,
                            });
                          } else {
                            calculateR('cosineSimilarityWithinPublic', {
                              profileType: withinProfileType,
                              matrixSize: withinMatrixSize,
                              study: study,
                              cancerType: cancerType,
                              experimentalStrategy: pExperimentalStrategy,
                            });
                          }
                        }}
                      >
                        Calculate
                      </Button>
                    </Col>
                  </Row>

                  <div id="withinPlot">
                    <div style={{ display: withinErr ? 'block' : 'none' }}>
                      <p>
                        An error has occured. Check the debug section for more
                        info.
                      </p>
                    </div>
                    <div style={{ display: withinPlotURL ? 'block' : 'none' }}>
                      <div className="d-flex">
                        <a
                          className="px-2 py-1"
                          href={withinPlotURL}
                          download={withinPlotURL.split('/').slice(-1)[0]}
                        >
                          Download Plot
                        </a>
                        <span className="ml-auto">
                          <Button
                            className="px-2 py-1"
                            variant="link"
                            onClick={() => downloadResults(withinTxtPath)}
                          >
                            Download Results
                          </Button>
                        </span>
                      </div>
                      <div className="p-2 border rounded">
                        <Row>
                          <Col>
                            <img
                              className="w-100 my-4 h-500"
                              src={withinPlotURL}
                            ></img>
                          </Col>
                        </Row>
                      </div>
                    </div>
                  </div>
                </div>
              </Form>
            </Body>
          </Collapse>
        </Card>
      </Accordion>

      <Accordion defaultActiveKey="1">
        <Card>
          <Toggle
            className="font-weight-bold"
            as={Header}
            eventKey="1"
            onClick={() =>
              dispatchCosineSimilarity({
                displayRefSig: !displayRefSig,
              })
            }
          >
            {displayRefSig == true ? (
              <FontAwesomeIcon icon={faMinus} />
            ) : (
              <FontAwesomeIcon icon={faPlus} />
            )}{' '}
            Cosine Similarity to Reference Signatures
          </Toggle>
          <Collapse eventKey="1">
            <Body>
              <Form className="my-2">
                <LoadingOverlay active={refSubmitOverlay} />
                <div>
                  <Row className="justify-content-center">
                    <Col sm="5">
                      <Group controlId="refProfileType">
                        <Label>Profile Type</Label>
                        <Select
                          options={profileOptions}
                          value={[refProfileType]}
                          onChange={(refProfileType) => {
                            dispatchCosineSimilarity({
                              refProfileType: refProfileType,
                            });
                            getSignatureSet(refProfileType);
                          }}
                          getOptionLabel={(option) => option}
                          getOptionValue={(option) => option}
                          {...selectFix}
                        />
                      </Group>
                    </Col>
                    <Col sm="5">
                      <Group controlId="refSignatureSet">
                        <Label>Reference Signature Set</Label>
                        <Select
                          options={refSignatureSetOptions}
                          value={[refSignatureSet]}
                          onChange={(refSignatureSet) => {
                            dispatchCosineSimilarity({
                              refSignatureSet: refSignatureSet,
                            });
                          }}
                          getOptionLabel={(option) => option}
                          getOptionValue={(option) => option}
                          {...selectFix}
                        />
                      </Group>
                    </Col>
                    <Col sm="2" className="m-auto">
                      <Button
                        variant="primary"
                        onClick={() => {
                          if (source == 'user') {
                            calculateR('cosineSimilarityRefSig', {
                              profileType: refProfileType,
                              signatureSet: refSignatureSet,
                              matrixList: JSON.stringify(
                                matrixList.filter(
                                  (matrix) =>
                                    matrix.Profile_Type == refProfileType
                                )
                              ),
                            });
                          } else {
                            calculateR('cosineSimilarityRefSigPublic', {
                              profileType: refProfileType,
                              signatureSet: refSignatureSet,
                              study: study,
                              cancerType: cancerType,
                              experimentalStrategy: pExperimentalStrategy,
                            });
                          }
                        }}
                      >
                        Calculate
                      </Button>
                    </Col>
                  </Row>

                  <div id="refPlot">
                    <div style={{ display: refErr ? 'block' : 'none' }}>
                      <p>
                        An error has occured. Check the debug section for more
                        info.
                      </p>
                    </div>
                    <div style={{ display: refPlotURL ? 'block' : 'none' }}>
                      <div className="d-flex">
                        <a
                          className="px-2 py-1"
                          href={refPlotURL}
                          download={refPlotURL.split('/').slice(-1)[0]}
                        >
                          Download Plot
                        </a>
                        <span className="ml-auto">
                          <Button
                            className="px-2 py-1"
                            variant="link"
                            onClick={() => downloadResults(refTxtPath)}
                          >
                            Download Results
                          </Button>
                        </span>
                      </div>
                      <div className="p-2 border rounded">
                        <Row>
                          <Col>
                            {refErr == true && (
                              <p>
                                An error has occured. Check the debug section
                                for more info.
                              </p>
                            )}
                            <img
                              className="w-100 my-4 h-500"
                              src={refPlotURL}
                            ></img>
                          </Col>
                        </Row>
                      </div>
                    </div>
                  </div>
                </div>
              </Form>
            </Body>
          </Collapse>
        </Card>
      </Accordion>

      <Button
        variant="link"
        className="p-0 mt-5"
        onClick={() =>
          dispatchCosineSimilarity({
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
