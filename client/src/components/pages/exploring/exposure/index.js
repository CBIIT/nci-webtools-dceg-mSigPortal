import React, { useState, useEffect } from 'react';
import { Form, Row, Col, Card, Button, Nav } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import Tumor from './tumor';
import Separated from './separated';
import Across from './across';
import Association from './association';
import Decomposition from './decomposition';
import Landscape from './landscape';
import Prevalence from './prevalence';
import {
  getInitialState,
  dispatchError,
  dispatchExploring,
  dispatchExpExposure,
  dispatchExpTumor,
  dispatchExpAcross,
  dispatchExpAssociation,
  dispatchExpDecomposition,
  dispatchExpLandscape,
  dispatchExpPrevalence,
  dispatchExpSeparated,
} from '../../../../services/store';
import Select from '../../../controls/select/select';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '../../../controls/sidebar-container/sidebar-container';
import Accordions from '../../../controls/accordions/accordions';

const { Header, Body } = Card;
const { Group, Label, Check } = Form;

export default function ExposureExploring({ populateControls }) {
  const { exposureSignature, exposureCancer, signatureNames } = useSelector(
    (state) => state.exploring
  );
  const { loading: loadingAcross } = useSelector((state) => state.expAcross);
  const { loading: loadingAssociation } = useSelector(
    (state) => state.expAssociation
  );
  const { loading: loadingLandscape } = useSelector(
    (state) => state.expLandscape
  );
  const { loading: loadingPrevalence } = useSelector(
    (state) => state.expPrevalence
  );

  const {
    study,
    studyOptions,
    strategy,
    strategyOptions,
    cancer,
    cancerOptions,
    refSigData,
    refSignatureSet,
    refSignatureSetOptions,
    signatureNameOptions,
    userNameOptions,
    genome,
    genomeOptions,
    exposureFile,
    matrixFile,
    signatureFile,
    usePublicSignature,
    source,
    loading,
    projectID,
    openSidebar,
  } = useSelector((state) => state.expExposure);
  const acrossArgs = useSelector((state) => state.expAcross);
  const associationArgs = useSelector((state) => state.expAssociation);
  const landscapeArgs = useSelector((state) => state.expLandscape);
  const prevalenceArgs = useSelector((state) => state.expPrevalence);

  const [exposureFileObj, setExposure] = useState(new File([], ''));
  const [matrixFileObj, setMatrix] = useState(new File([], ''));
  const [signatureFileObj, setSignature] = useState(new File([], ''));
  const [variableFileObj, setVariable] = useState(new File([], ''));

  const [exposureValidity, setExposureValidity] = useState(false);
  const [matrixValidity, setMatrixValidity] = useState(false);
  const [signatureValidity, setSignatureValidity] = useState(false);

  useEffect(() => {
    dispatchExploring({ displayTab: 'exposure' });
  }, []);

  function submitR(fn, args, id = projectID) {
    return fetch(`api/exploringR`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        fn: fn,
        args: args,
        projectID: id,
      }),
    }).then((res) => res.json());
  }

  async function calculateAcross() {
    if (source == 'user' && !projectID) {
      dispatchError('Missing Required Files');
    } else {
      dispatchExpAcross({
        loading: true,
        err: false,
        debugR: '',
      });

      if (source == 'user') {
        if (!projectID) {
          try {
            const id = handleUpload();
            await handleCalculate('across', id);
          } catch (error) {
            dispatchError(error);
          }
        }
      } else {
        await handleCalculate('across');
      }

      dispatchExpAcross({ loading: false });
    }
  }

  async function calculateAssociation() {
    if (source == 'user' && !projectID) {
      dispatchError('Missing Required Files');
    } else {
      dispatchExpAssociation({
        loading: true,
        err: false,
        debugR: '',
      });

      if (source == 'user') {
        if (!projectID) {
          try {
            const id = await handleUpload();
            await handleCalculate('association', id);
          } catch (error) {
            dispatchError(error);
          }
        }
      } else {
        await handleCalculate('association');
      }

      dispatchExpAssociation({ loading: false });
    }
  }

  async function calculateLandscape() {
    if (source == 'user' && !projectID) {
      dispatchError('Missing Required Files');
    } else {
      dispatchExpLandscape({
        loading: true,
        err: false,
        debugR: '',
      });

      if (source == 'user') {
        if (!projectID) {
          try {
            const id = await handleUpload();
            await handleCalculate('landscape', id);
          } catch (error) {
            dispatchError(error);
          }
        }
      } else if (variableFileObj.size) {
        try {
          const id = await uploadVariable();
          await handleCalculate('landscape', id);
        } catch (error) {
          dispatchError(error);
        }
      } else {
        await handleCalculate('landscape');
      }

      dispatchExpLandscape({ loading: false });
    }
  }

  async function calculatePrevalence() {
    if (source == 'user' && !projectID) {
      dispatchError('Missing Required Files');
    } else {
      dispatchExpPrevalence({
        loading: true,
        err: false,
        debugR: '',
      });

      if (source == 'user') {
        if (!projectID) {
          try {
            const id = await handleUpload();
            await handleCalculate('prevalence', id);
          } catch (error) {
            dispatchError(error);
          }
        } else {
          await handleCalculate('prevalence', projectID);
        }
      } else {
        await handleCalculate('prevalence');
      }

      dispatchExpPrevalence({ loading: false });
    }
  }

  async function calculateAll() {
    try {
      dispatchExpExposure({ loading: true });
      if (source == 'user') {
        const { projectID, exposureData } = await handleUpload();

        // get signature name options, ignore sample key
        const nameOptions = Object.keys(exposureData[0]).filter(
          (key) => key != 'Samples'
        );

        dispatchExpAcross({ signatureName: nameOptions[0] });
        dispatchExpAssociation({
          signatureName1: nameOptions[0],
          signatureName2: nameOptions[1],
        });
        dispatchExpExposure({
          projectID: projectID,
          userNameOptions: nameOptions,
        });
        try {
          await handleCalculate('all', projectID);
        } catch (error) {
          dispatchError(error);
        }
      } else if (variableFileObj.size) {
        try {
          const id = await uploadVariable();
          await handleCalculate('all', id);
        } catch (error) {
          dispatchError(error);
        }
      } else {
        await handleCalculate('all');
      }
    } catch (err) {
      dispatchError(err);
    } finally {
      dispatchExpExposure({ loading: false });
    }
  }

  async function handleCalculate(fn = 'all', id = projectID) {
    let rFn = 'exposurePublic';
    let args = {
      fn: fn,
      common: JSON.stringify({
        study: study,
        strategy: strategy,
        refSignatureSet: refSignatureSet,
        cancerType: cancer,
        genome: genome,
      }),
    };
    if (!loadingAcross && (fn == 'all' || fn == 'across')) {
      args.across = JSON.stringify({
        signatureName: acrossArgs.signatureName,
      });
    }
    if (!loadingAssociation && (fn == 'all' || fn == 'association')) {
      args.association = JSON.stringify({
        useCancerType: associationArgs.toggleCancer,
        both: associationArgs.both,
        signatureName1: associationArgs.signatureName1,
        signatureName2: associationArgs.signatureName2,
      });
    }
    if (!loadingLandscape && (fn == 'all' || fn == 'landscape')) {
      args.landscape = JSON.stringify({
        variableFile: landscapeArgs.variableFile,
      });
    }
    if (!loadingPrevalence && (fn == 'all' || fn == 'prevalence')) {
      args.prevalence = JSON.stringify({
        mutation: parseFloat(prevalenceArgs.mutation) || 100,
      });
    }
    if (source == 'user') {
      rFn = 'exposureUser';
      args.files = JSON.stringify({
        exposureFile: exposureFile,
        matrixFile: matrixFile,
        signatureFile: signatureFile,
      });
    }
    const { debugR, output, errors } = await submitR(rFn, args, id);

    if (output) {
      if (output.tumorPath)
        dispatchExpTumor({
          plotPath: output.tumorPath,
          debugR: debugR,
          err: false,
        });
      else if (fn == 'all')
        dispatchExpTumor({
          err: errors.tumorError,
          debugR: debugR,
          plotPath: '',
        });

      if (output.burdenSeparatedPath)
        dispatchExpSeparated({
          plotPath: output.burdenSeparatedPath,
          debugR: debugR,
          err: false,
        });
      else if (fn == 'all')
        dispatchExpSeparated({
          err: errors.burdenSeparatedError,
          debugR: debugR,
          plotPath: '',
        });

      if (output.burdenAcrossPath)
        dispatchExpAcross({
          plotPath: output.burdenAcrossPath,
          debugR: debugR,
          err: false,
        });
      else if (fn == 'all' || fn == 'across')
        dispatchExpAcross({
          err: errors.burdenAcrossError,
          debugR: debugR,
          plotPath: '',
        });

      if (output.associationPath)
        dispatchExpAssociation({
          plotPath: output.associationPath,
          debugR: debugR,
          err: false,
        });
      else if (fn == 'all' || fn == 'association')
        dispatchExpAssociation({
          err: errors.associationError,
          debugR: debugR,
          plotPath: '',
        });

      if (output.decompositionPath)
        dispatchExpDecomposition({
          plotPath: output.decompositionPath,
          txtPath: output.decompositionData,
          debugR: debugR,
          err: false,
        });
      else if (fn == 'all')
        dispatchExpDecomposition({
          err: errors.decompositionError,
          debugR: debugR,
          plotPath: '',
        });

      if (output.landscapePath)
        dispatchExpLandscape({
          plotPath: output.landscapePath,
          debugR: debugR,
          err: false,
        });
      else if (fn == 'all' || fn == 'landscape')
        dispatchExpLandscape({
          err: errors.landscapeError,
          debugR: debugR,
          plotPath: '',
        });

      if (output.prevalencePath)
        dispatchExpPrevalence({
          plotPath: output.prevalencePath,
          debugR: debugR,
          err: false,
        });
      else if (fn == 'all' || fn == 'prevalence')
        dispatchExpPrevalence({
          err: errors.prevalenceError,
          debugR: debugR,
          plotPath: '',
        });
    } else {
      dispatchError(debugR);
    }
  }

  function handleStudy(study) {
    const strategyOptions = [
      ...new Set(
        exposureSignature
          .filter((data) => data.Study == study)
          .map((data) => data.Dataset)
      ),
    ];
    const strategy = strategyOptions[0];

    const refSignatureSetOptions = [
      ...new Set(
        exposureSignature
          .filter((row) => row.Study == study && row.Dataset == strategy)
          .map((row) => row.Signature_set_name)
      ),
    ];
    const refSignatureSet = refSignatureSetOptions[0];

    const cancerOptions = [
      ...new Set(
        exposureCancer
          .filter((data) => data.Study == study && data.Dataset == strategy)
          .map((data) => data.Cancer_Type)
      ),
    ];

    handleSet(refSignatureSet);

    dispatchExpExposure({
      study: study,
      strategy: strategy,
      strategyOptions: strategyOptions,
      cancer: cancerOptions[0],
      cancerOptions: cancerOptions,
      refSignatureSetOptions: refSignatureSetOptions,
      refSignatureSet: refSignatureSet,
    });
  }

  function handleStrategy(strategy) {
    const refSignatureSetOptions = [
      ...new Set(
        exposureSignature
          .filter((row) => row.Study == study && row.Dataset == strategy)
          .map((row) => row.Signature_set_name)
      ),
    ];
    const refSignatureSet = refSignatureSetOptions[0];

    const cancerOptions = [
      ...new Set(
        exposureCancer
          .filter((data) => data.Study == study && data.Dataset == strategy)
          .map((data) => data.Cancer_Type)
      ),
    ];

    handleSet(refSignatureSet);

    dispatchExpExposure({
      strategy: strategy,
      cancer: cancerOptions[0],
      cancerOptions: cancerOptions,
      refSignatureSetOptions: refSignatureSetOptions,
      refSignatureSet: refSignatureSet,
    });
  }

  function handleSet(set) {
    const signatureNameOptions = [
      ...new Set(
        signatureNames
          .filter((row) => row.Signature_set_name == set)
          .map((row) => row.Signature_name)
      ),
    ];

    dispatchExpExposure({
      refSignatureSet: set,
      signatureNameOptions: signatureNameOptions,
    });
  }

  async function handleUpload() {
    return new Promise(async (resolve, reject) => {
      if (
        exposureFileObj.size &&
        matrixFileObj.size &&
        ((!usePublicSignature && signatureFileObj) ||
          (usePublicSignature && refSignatureSet))
      ) {
        try {
          const data = new FormData();
          data.append('exposureFile', exposureFileObj);
          data.append('matrixFile', matrixFileObj);
          if (!usePublicSignature)
            data.append('signatureFile', signatureFileObj);
          if (variableFileObj.size)
            data.append('variableFile', variableFileObj);
          let response = await fetch(`api/upload`, {
            method: 'POST',
            body: data,
          });

          if (!response.ok) {
            const { msg, error } = await response.json();
            const message = `<div>
            <p>${msg}</p>
          ${error ? `<p>${error}</p>` : ''} 
          </div>`;
            dispatchError(message);
            reject(error);
          } else {
            const { projectID, exposurePath } = await response.json();

            const exposureData = await (
              await fetch('api/getSignaturesUser', {
                method: 'POST',
                headers: {
                  Accept: 'application/json',
                  'Content-Type': 'application/json',
                },
                body: JSON.stringify({
                  path: exposurePath,
                }),
              })
            ).json();
            resolve({ projectID, exposureData });
          }
        } catch (err) {
          dispatchError(err);
          reject(err);
        }
      } else {
        reject('Missing required files');
      }
    });
  }

  function handleReset() {
    const initialState = getInitialState();

    source == 'public'
      ? dispatchExpExposure(initialState.expExposure)
      : dispatchExpExposure({ ...initialState.expExposure, source: 'user' });
    dispatchExpTumor(initialState.expTumor);
    dispatchExpSeparated(initialState.expSeparated);
    dispatchExpAcross(initialState.expAcross);
    dispatchExpDecomposition(initialState.expDecomposition);
    dispatchExpAssociation(initialState.expAssociation);
    dispatchExpLandscape(initialState.expLandscape);
    dispatchExpPrevalence(initialState.expPrevalence);
    populateControls();
  }

  // when using public data and only need to upload a variable data file
  async function uploadVariable() {
    return new Promise(async (resolve, reject) => {
      try {
        const data = new FormData();
        data.append('variableFile', variableFileObj);

        let response = await fetch(`api/upload`, {
          method: 'POST',
          body: data,
        });

        if (!response.ok) {
          const { msg, error } = await response.json();
          const message = `<div>
                            <p>${msg}</p>
                            ${error ? `<p>${error}</p>` : ''} 
                          </div>`;
          dispatchError(message);
          reject(error);
        } else {
          const { projectID } = await response.json();
          resolve(projectID);
        }
      } catch (err) {
        dispatchError(err);
        reject(err);
      }
    });
  }

  function handleVariable(file) {
    setVariable(file);
    dispatchExpLandscape({ variableFile: file.name });
  }

  const sections = [
    {
      component: <Tumor />,
      id: 'tumor',
      title: 'Tumor Mutational Burden',
    },
    {
      component: <Separated />,
      id: 'separated',
      title: 'Tumor Mutational Burden Separated by Signatures',
    },
    {
      component: <Across calculateAcross={calculateAcross} />,
      id: 'across',
      title: 'Mutational Signature Burden Across Cancer Types',
    },
    {
      component: <Decomposition />,
      id: 'decomposition',
      title: 'Evaluating the Performance of Mutational Signature Decomposition',
    },
    {
      component: <Association calculateAssociation={calculateAssociation} />,
      id: 'association',
      title: 'Mutational Signature Association',
    },
    {
      component: (
        <Landscape
          calculateLandscape={calculateLandscape}
          handleVariable={handleVariable}
        />
      ),
      id: 'landscape',
      title: 'Landscape of Mutational Signature Activity',
    },
    {
      component: <Prevalence calculatePrevalence={calculatePrevalence} />,
      id: 'prevalence',
      title: 'Prevalence of Mutational Signature',
    },
  ];

  return (
    <div className="position-relative">
      <SidebarContainer
        collapsed={!openSidebar}
        onCollapsed={(e) => dispatchExpExposure({ openSidebar: !e })}
      >
        <SidebarPanel>
          <div className="p-3 bg-white border rounded">
            <Row>
              <Col sm="auto">
                <Group>
                  <Label className="mr-auto">
                    <h3 className="mb-2">Data Source</h3>
                  </Label>
                  <Check inline id="radioPublic" className="ml-4">
                    <Check.Input
                      disabled={loading}
                      type="radio"
                      value="public"
                      checked={source == 'public'}
                      onChange={(e) =>
                        dispatchExpExposure({ source: 'public' })
                      }
                    />
                    <Check.Label className="font-weight-normal">
                      Public
                    </Check.Label>
                  </Check>
                  <Check inline id="radioUser">
                    <Check.Input
                      disabled={loading}
                      type="radio"
                      value="user"
                      checked={source == 'user'}
                      onChange={(e) => dispatchExpExposure({ source: 'user' })}
                    />
                    <Check.Label className="font-weight-normal">
                      User
                    </Check.Label>
                  </Check>
                </Group>
              </Col>
            </Row>
            {source == 'public' ? (
              <div>
                <Row>
                  <Col>
                    <Group>
                      <Select
                        disabled={loading}
                        id="expStudyPublic"
                        label="Study"
                        value={study}
                        options={studyOptions}
                        onChange={handleStudy}
                      />
                    </Group>
                  </Col>
                </Row>
                <Row>
                  <Col>
                    <Group>
                      <Select
                        disabled={loading}
                        id="tumorStrategy"
                        label="Experimental Strategy"
                        value={strategy}
                        options={strategyOptions}
                        onChange={handleStrategy}
                      />
                    </Group>
                  </Col>
                </Row>
                <Row>
                  <Col>
                    <Group>
                      <Select
                        disabled={loading}
                        id="expSetPublic"
                        label="Reference Signature Set"
                        value={refSignatureSet}
                        options={refSignatureSetOptions}
                        onChange={handleSet}
                      />
                    </Group>
                  </Col>
                </Row>
                <Row>
                  <Col>
                    <Group>
                      <Select
                        className="mb-4"
                        disabled={loading}
                        id="prevalenceCancerType"
                        label="Cancer Type"
                        value={cancer}
                        options={cancerOptions}
                        onChange={(cancer) =>
                          dispatchExpExposure({
                            cancer: cancer,
                          })
                        }
                      />
                    </Group>
                  </Col>
                </Row>
                <Row>
                  <Col sm="12" md="6">
                    <Button
                      disabled={loading}
                      className="w-100"
                      variant="secondary"
                      onClick={() => handleReset()}
                    >
                      Reset
                    </Button>
                  </Col>
                  <Col sm="12" md="6">
                    <Button
                      disabled={loading}
                      className="w-100"
                      variant="primary"
                      onClick={() => calculateAll()}
                    >
                      Calculate
                    </Button>
                  </Col>
                </Row>
              </div>
            ) : (
              <div>
                <Row>
                  <Col>
                    <Group>
                      <Label>Upload Exposure File</Label>
                      <Form.File
                        disabled={loading}
                        id="variableData"
                        label={exposureFileObj.name || 'Exposure File'}
                        accept=".txt"
                        onChange={(e) => {
                          setExposure(e.target.files[0]);
                          dispatchExpExposure({
                            exposureFile: e.target.files[0].name,
                          });
                        }}
                        custom
                      />
                      {exposureValidity && (
                        <span className="text-danger">
                          Exposure File Required
                        </span>
                      )}
                    </Group>
                  </Col>
                </Row>
                <Row>
                  <Col>
                    <Group>
                      <Label>Upload Matrix File</Label>
                      <Form.File
                        disabled={loading}
                        id="variableData"
                        label={matrixFileObj.name || 'Matrix File'}
                        accept=".txt"
                        onChange={(e) => {
                          setMatrix(e.target.files[0]);
                          dispatchExpExposure({
                            matrixFile: e.target.files[0].name,
                          });
                        }}
                        custom
                      />
                      {matrixValidity && (
                        <span className="text-danger">
                          Matrix File Required
                        </span>
                      )}
                    </Group>
                  </Col>
                </Row>
                <Row>
                  <Col>
                    <Group controlId="toggleSignatureSource" className="d-flex">
                      <Label className="mr-4">Use Public Signature Data</Label>
                      <Check inline id="toggleSignatureSource">
                        <Check.Input
                          disabled={loading}
                          type="checkbox"
                          value={usePublicSignature}
                          checked={usePublicSignature}
                          onChange={() =>
                            dispatchExpExposure({
                              usePublicSignature: !usePublicSignature,
                            })
                          }
                        />
                      </Check>
                    </Group>
                  </Col>
                </Row>

                {usePublicSignature ? (
                  <div>
                    <Row>
                      <Col>
                        <Group>
                          <Select
                            disabled={loading}
                            id="expStudyUser"
                            label="Study"
                            value={study}
                            options={studyOptions}
                            onChange={handleStudy}
                          />
                        </Group>
                      </Col>
                    </Row>
                    <Row>
                      <Col>
                        <Group>
                          <Select
                            disabled={loading}
                            id="exposureSignatureSet"
                            label="Reference Signature Set"
                            value={refSignatureSet}
                            options={refSignatureSetOptions}
                            onChange={handleSet}
                          />{' '}
                        </Group>
                      </Col>
                    </Row>
                  </div>
                ) : (
                  <Row>
                    <Col>
                      <Group>
                        <Label>Upload Signature Data</Label>
                        <Form.File
                          disabled={loading}
                          id="variableData"
                          label={signatureFileObj.name || 'Signature File'}
                          accept=".txt"
                          onChange={(e) => {
                            setSignature(e.target.files[0]);
                            dispatchExpExposure({
                              signatureFile: e.target.files[0].name,
                            });
                          }}
                          custom
                        />
                        {signatureValidity && (
                          <span className="text-danger">
                            Signature File Required
                          </span>
                        )}
                      </Group>
                    </Col>
                  </Row>
                )}
                <Row>
                  <Col>
                    <Group>
                      <Select
                        disabled={loading}
                        id="exposureGenome"
                        label="Genome"
                        value={genome}
                        options={genomeOptions}
                        onChange={(genome) =>
                          dispatchExpExposure({ genome: genome })
                        }
                      />
                    </Group>
                  </Col>
                </Row>
                <Row>
                  <Col sm="6">
                    <Button
                      disabled={loading}
                      className="w-100"
                      variant="secondary"
                      onClick={() => handleReset()}
                    >
                      Reset
                    </Button>
                  </Col>
                  <Col sm="6">
                    <Button
                      disabled={loading}
                      className="w-100"
                      variant="primary"
                      onClick={() => calculateAll()}
                    >
                      Calculate
                    </Button>
                  </Col>
                </Row>
              </div>
            )}
          </div>
        </SidebarPanel>
        <MainPanel>
          <Card>
            <Header>
              <Nav variant="pills" defaultActiveKey="#exploring/exposure">
                <Nav.Item>
                  <Nav.Link href="#exploring/exposure">
                    Exposure Exploring
                  </Nav.Link>
                </Nav.Item>
              </Nav>
            </Header>
            <Body>
              <LoadingOverlay active={loading} />
              <Accordions components={sections} />
            </Body>
          </Card>
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
