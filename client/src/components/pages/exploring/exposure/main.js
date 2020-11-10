import React, { useState } from 'react';
import { Form, Row, Col, Accordion, Card, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import Tumor from './tumor';
import Activity from './activity';
import Association from './association';
import Decomposition from './decomposition';
import Landscape from './landscape';
import Prevalence from './prevalence';
import {
  dispatchError,
  dispatchExploring,
  dispatchExpExposure,
  dispatchExpTumor,
  dispatchExpActivity,
  dispatchExpAssociation,
  dispatchExpDecomposition,
  dispatchExpLandscape,
  dispatchExpPrevalence,
} from '../../../../services/store';
import Select from '../../../controls/select/select';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '../../../controls/sidebar-container/sidebar-container';

const { Header, Body } = Card;
const { Toggle, Collapse } = Accordion;
const { Group, Label, Check, Control } = Form;

export default function ExposureExploring() {
  const rootURL = window.location.pathname;
  const { displayTab, exposureAccordion, publicDataOptions } = useSelector(
    (state) => state.exploring
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
  const activityArgs = useSelector((state) => state.expActivity);
  const associationArgs = useSelector((state) => state.expAssociation);
  const landscapeArgs = useSelector((state) => state.expLandscape);
  const prevalenceArgs = useSelector((state) => state.expPrevalence);

  const [exposureFileObj, setExposure] = useState(new File([], ''));
  const [matrixFileObj, setMatrix] = useState(new File([], ''));
  const [signatureFileObj, setSignature] = useState(new File([], ''));
  const [exposureValidity, setExposureValidity] = useState(false);
  const [matrixValidity, setMatrixValidity] = useState(false);
  const [signatureValidity, setSignatureValidity] = useState(false);

  function submitR(fn, args, id = projectID) {
    return fetch(`${rootURL}exploringR`, {
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

  async function calculateActivity() {
    dispatchExpActivity({
      loading: true,
      err: false,
      debugR: '',
    });

    if (source == 'user') {
      if (!projectID) {
        try {
          const id = await handleUpload();
          await handleCalculate('activity', id);
        } catch (error) {
          dispatchError(error);
        }
      }
    } else {
      await handleCalculate('activity');
    }

    dispatchExpActivity({ loading: false });
  }

  async function calculateAssociation() {
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

  async function calculateLandscape() {
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
    } else {
      await handleCalculate('landscape');
    }

    dispatchExpLandscape({ loading: false });
  }

  async function calculatePrevalence() {
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

  async function calculateAll() {
    dispatchExpExposure({ loading: true });

    if (source == 'user') {
      if (!projectID) {
        try {
          const id = await handleUpload();
          await handleCalculate('all', id);
        } catch (error) {
          dispatchError(error);
        }
      }
    } else {
      await handleCalculate('all');
    }

    dispatchExpExposure({ loading: false });
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
    if (fn == 'all' || fn == 'activity') {
      args.activity = JSON.stringify({
        signatureName: activityArgs.signatureName,
      });
    }
    if (fn == 'all' || fn == 'association') {
      args.association = JSON.stringify({
        useCancerType: associationArgs.toggleCancer,
        both: associationArgs.both,
        signatureName1: associationArgs.signatureName1,
        signatureName2: associationArgs.signatureName2,
      });
    }
    if (fn == 'all' || fn == 'landscape') {
      args.landscape = JSON.stringify({
        varDataPath: landscapeArgs.varDataPath,
      });
    }
    if (fn == 'all' || fn == 'prevalence') {
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
    const { debugR, output } = await submitR(rFn, args, id);

    if (output) {
      if (output.tumorPath)
        dispatchExpTumor({
          plotPath: output.tumorPath,
          debugR: debugR,
          err: false,
        });
      else if (fn == 'all') dispatchExpTumor({ err: true, debugR: debugR });

      if (output.activityPath)
        dispatchExpActivity({
          plotPath: output.activityPath,
          debugR: debugR,
          err: false,
        });
      else if (fn == 'all' || fn == 'activity')
        dispatchExpActivity({ err: true, debugR: debugR });

      if (output.associationPath)
        dispatchExpAssociation({
          plotPath: output.associationPath,
          debugR: debugR,
          err: false,
        });
      else if (fn == 'all' || fn == 'association')
        dispatchExpAssociation({ err: true, debugR: debugR });

      if (output.decompositionPath)
        dispatchExpDecomposition({
          plotPath: output.decompositionPath,
          txtPath: output.decompositionData,
          debugR: debugR,
          err: false,
        });
      else if (fn == 'all')
        dispatchExpDecomposition({ err: true, debugR: debugR });

      if (output.landscapePath)
        dispatchExpLandscape({
          plotPath: output.landscapePath,
          debugR: debugR,
          err: false,
        });
      else if (fn == 'all' || fn == 'landscape')
        dispatchExpLandscape({ err: true, debugR: debugR });

      if (output.prevalencePath)
        dispatchExpPrevalence({
          plotPath: output.prevalencePath,
          debugR: debugR,
          err: false,
        });
      else if (fn == 'all' || fn == 'prevalence')
        dispatchExpPrevalence({ err: true, debugR: debugR });
    } else {
      dispatchError(debugR);
    }
  }

  function handleStudy(study) {
    const strategyOptions = [
      ...new Set(
        publicDataOptions
          .filter((data) => data.Study == study)
          .map((data) => data.Dataset)
      ),
    ];

    const cancerOptions = [
      ...new Set(
        publicDataOptions
          .filter((data) => data.Study == study)
          .map((data) => data.Cancer_Type)
      ),
    ];

    dispatchExpExposure({
      study: study,
      strategy: strategyOptions[0],
      strategyOptions: strategyOptions,
      cancer: cancerOptions[0],
      cancerOptions: cancerOptions,
    });
  }

  function handleSet(set) {
    const signatureNameOptions = [
      ...new Set(
        refSigData
          .filter((row) => row.Signature_set_name == set)
          .map((row) => row.Signature_name)
      ),
    ];
    dispatchExpExposure({
      refSignatureSet: set,
      signatureNameOptions: signatureNameOptions,
    });
  }

  function handleUpload() {
    return new Promise(async (resolve, reject) => {
      if (
        exposureFileObj.size &&
        matrixFileObj.size &&
        ((!usePublicSignature && signatureFileObj) ||
          (usePublicSignature && refSignatureSet))
      ) {
        try {
          const data = new FormData();
          data.append('inputFile', exposureFileObj);
          data.append('inputFile', matrixFileObj);
          if (!usePublicSignature) data.append('inputFile', signatureFileObj);
          let response = await fetch(`${rootURL}upload`, {
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
      } else {
        reject('Missing required files');
      }
    });
  }

  const sections = [
    {
      component: <Tumor />,
      id: 'tumor',
      title: 'Tumor Mutational Burden',
    },
    {
      component: <Activity calculateActivity={calculateActivity} />,
      id: 'activity',
      title: 'Mutational Signature Activity',
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
      component: <Landscape calculateLandscape={calculateLandscape} />,
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
            <LoadingOverlay active={loading} />
            <Row>
              <Col sm="auto">
                <Group>
                  <Label className="mr-auto">
                    <h3 className="mb-2">Data Source</h3>
                  </Label>
                  <Check inline id="radioPublic" className="ml-4">
                    <Check.Input
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
                        id="tumorStudy"
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
                        id="tumorStrategy"
                        label="Experimental Strategy"
                        value={strategy}
                        options={strategyOptions}
                        onChange={(strategy) =>
                          dispatchExpExposure({ strategy: strategy })
                        }
                      />
                    </Group>
                  </Col>
                </Row>
                <Row>
                  <Col>
                    <Group>
                      <Select
                        id="prevalenceCancerType"
                        label="Cancer Type"
                        value={cancer}
                        options={cancerOptions}
                        onChange={(cancer) =>
                          dispatchExpPrevalence({
                            cancer: cancer,
                          })
                        }
                      />
                    </Group>
                  </Col>
                </Row>
                <Row>
                  <Col>
                    <Group>
                      <Select
                        id="tumorSet"
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
                      <Button variant="primary" onClick={() => calculateAll()}>
                        Calculate
                      </Button>
                    </Group>
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

                <Row>
                  <Col>
                    <Group>
                      {usePublicSignature ? (
                        <Select
                          id="exposureSignatureSet"
                          label="Reference Signature Set"
                          value={refSignatureSet}
                          options={refSignatureSetOptions}
                          onChange={handleSet}
                        />
                      ) : (
                        <div>
                          <Label>Upload Signature Data</Label>
                          <Form.File
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
                        </div>
                      )}
                    </Group>
                  </Col>
                </Row>
                <Row>
                  <Col>
                    <Group>
                      <Select
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
                  <Col>
                    <Group>
                      <Button variant="primary" onClick={() => calculateAll()}>
                        Calculate
                      </Button>
                    </Group>
                  </Col>
                </Row>
              </div>
            )}
          </div>
        </SidebarPanel>
        <MainPanel>
          {sections.map(({ component, id, title }) => {
            return (
              <Accordion activeKey={exposureAccordion[id]} key={id}>
                <Card>
                  <Toggle
                    className="font-weight-bold"
                    as={Header}
                    eventKey={exposureAccordion[id]}
                    onClick={() =>
                      dispatchExploring({
                        exposureAccordion: {
                          ...exposureAccordion,
                          [id]: !exposureAccordion[id],
                        },
                      })
                    }
                  >
                    {exposureAccordion[id] == true ? (
                      <FontAwesomeIcon icon={faMinus} />
                    ) : (
                      <FontAwesomeIcon icon={faPlus} />
                    )}{' '}
                    {title}
                  </Toggle>
                  <Collapse eventKey={exposureAccordion[id]}>
                    <Body>{component}</Body>
                  </Collapse>
                </Card>
              </Accordion>
            );
          })}
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
