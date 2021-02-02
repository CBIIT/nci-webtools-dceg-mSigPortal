import React, { useState, useEffect, useRef } from 'react';
import { Form, Row, Col, Tab, Button, Nav } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import TMB from './tmb';
import TmbSig from './tmbSignatures';
import MsBurden from './msBurden';
import MsAssociation from './msAssociation';
import MsDecomposition from './msDecomposition';
import MsLandscape from './msLandscape';
import MsPrevalence from './msPrevalence';
import {
  getInitialState,
  dispatchError,
  dispatchExploring,
  dispatchExpExposure,
  dispatchTMB,
  dispatchTmbSignatures,
  dispatchMsBurden,
  dispatchMsAssociation,
  dispatchMsDecomposition,
  dispatchMsLandscape,
  dispatchMsPrevalence,
} from '../../../../services/store';
import Select from '../../../controls/select/select';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '../../../controls/sidebar-container/sidebar-container';

const { Group, Label, Check } = Form;
const { Container, Content, Pane } = Tab;
const { Item, Link } = Nav;

export default function Exposure({ match, populateControls }) {
  const { exampleName } = match.params;
  const { exposureSignature, exposureCancer, signatureNames } = useSelector(
    (state) => state.exploring
  );
  const { loading: loadingMsBurden } = useSelector((state) => state.msBurden);
  const { loading: loadingMsAssociation } = useSelector(
    (state) => state.msAssociation
  );
  const { loading: loadingMsLandscape } = useSelector(
    (state) => state.msLandscape
  );
  const { loading: loadingMsPrevalence } = useSelector(
    (state) => state.msPrevalence
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
    display,
    loading,
    loadingMsg,
    projectID,
    openSidebar,
  } = useSelector((state) => state.expExposure);
  const burdenArgs = useSelector((state) => state.msBurden);
  const associationArgs = useSelector((state) => state.msAssociation);
  const landscapeArgs = useSelector((state) => state.msLandscape);
  const prevalenceArgs = useSelector((state) => state.msPrevalence);

  const [exposureFileObj, setExposure] = useState(new File([], ''));
  const [matrixFileObj, setMatrix] = useState(new File([], ''));
  const [signatureFileObj, setSignature] = useState(new File([], ''));
  const [variableFileObj, setVariable] = useState(new File([], ''));

  const [exposureValidity, setExposureValidity] = useState(false);
  const [matrixValidity, setMatrixValidity] = useState(false);
  const [signatureValidity, setSignatureValidity] = useState(false);

  // load example if available
  useEffect(() => {
    if (exampleName) loadExample(exampleName);
  }, [exampleName]);

  // set selected tab on component render
  useEffect(() => {
    dispatchExploring({ displayTab: 'exposure' });
  }, []);

  // get new signature name options filtered by cancer type on first mount
  function usePrevious(value) {
    const ref = useRef();
    useEffect(() => {
      ref.current = value;
    }, [value]);

    return ref.current;
  }
  const prevSigNameOptions = usePrevious(signatureNameOptions);
  useEffect(() => {
    if (
      associationArgs.toggleCancer &&
      study &&
      strategy &&
      refSignatureSet &&
      source == 'public' &&
      prevSigNameOptions
    )
      getSignatureNames();
  }, [study, strategy, refSignatureSet, cancer, associationArgs.toggleCancer]);

  async function loadExample(id) {
    dispatchExpExposure({
      loading: {
        active: true,
        content: 'Loading Example',
        showIndicator: true,
      },
    });
    try {
      const { projectID, state } = await (
        await fetch(`api/getExposureExample/${id}`)
      ).json();

      dispatchExpExposure({ ...state.expExposure, projectID: projectID });
      // rehydrate state if available
      if (state.tmb) dispatchTMB(state.tmb);
      if (state.msBurden) dispatchMsBurden(state.msBurden);
      if (state.msAssociation) dispatchMsAssociation(state.msAssociation);
      if (state.msDecomposition) dispatchMsDecomposition(state.msDecomposition);
      if (state.msLandscape) dispatchMsLandscape(state.msLandscape);
      if (state.msPrevalence) dispatchMsPrevalence(state.msPrevalence);
      if (state.tmbSignatures) dispatchTmbSignatures(state.tmbSignatures);
    } catch (error) {
      dispatchError(error);
    }
    dispatchExpExposure({
      loading: false,
    });
  }

  // get signature name options filtered by cancer type
  async function getSignatureNames() {
    dispatchExpExposure({
      loading: true,
      loadingMsg: 'Filtering Signature Names',
    });
    try {
      const { stdout, output } = await (
        await fetch(`api/getSignatureNames`, {
          method: 'POST',
          headers: {
            Accept: 'application/json',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            args: {
              study: study,
              strategy: strategy,
              refSignatureSet: refSignatureSet,
              cancerType: cancer,
            },
          }),
        })
      ).json();

      if (output.data.length)
        dispatchExpExposure({
          signatureNameOptions: output.data,
        });
      else dispatchError(stdout);
    } catch (err) {
      dispatchError(err);
    }
    dispatchExpExposure({ loading: false, loadingMsg: null });
  }

  async function submitR(fn, args, id = projectID) {
    try {
      const response = await fetch(`api/exploringR`, {
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
      });

      if (response.ok) {
        return await response.json();
      } else {
        dispatchError('R submit failed');
      }
    } catch (err) {
      dispatchError(err);
    }
  }

  async function calculateBurden() {
    if (source == 'user' && !projectID) {
      dispatchError('Missing Required Files');
    } else {
      dispatchMsBurden({
        loading: true,
        err: false,
        plotPath: '',
      });

      if (source == 'user') {
        if (!projectID) {
          try {
            const id = handleUpload();
            await handleCalculate('burden', id);
          } catch (error) {
            dispatchError(error);
          }
        }
      } else {
        await handleCalculate('burden');
      }

      dispatchMsBurden({ loading: false });
    }
  }

  async function calculateAssociation() {
    if (source == 'user' && !projectID) {
      dispatchError('Missing Required Files');
    } else {
      dispatchMsAssociation({
        loading: true,
        err: false,
        plotPath: '',
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

      dispatchMsAssociation({ loading: false });
    }
  }

  async function calculateLandscape() {
    if (source == 'user' && !projectID) {
      dispatchError('Missing Required Files');
    } else {
      dispatchMsLandscape({
        loading: true,
        err: false,
        plotPath: '',
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

      dispatchMsLandscape({ loading: false });
    }
  }

  async function calculatePrevalence() {
    if (source == 'user' && !projectID) {
      dispatchError('Missing Required Files');
    } else {
      dispatchMsPrevalence({
        loading: true,
        err: false,
        plotPath: '',
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

      dispatchMsPrevalence({ loading: false });
    }
  }

  async function calculateAll() {
    try {
      dispatchExpExposure({ loading: true });

      dispatchTMB({
        plotPath: '',
        err: '',
      });
      dispatchTmbSignatures({
        plotPath: '',
        err: '',
      });
      dispatchMsDecomposition({
        plotPath: '',
        txtPath: '',
        err: '',
      });
      dispatchMsBurden({
        plotPath: '',
        err: '',
      });
      dispatchMsAssociation({
        plotPath: '',
        err: '',
      });
      dispatchMsLandscape({
        plotPath: '',
        err: '',
      });

      dispatchMsPrevalence({
        plotPath: '',
        err: '',
      });

      if (source == 'user') {
        const { projectID, exposureData } = await handleUpload();

        // get signature name options, ignore sample key
        const nameOptions = Object.keys(exposureData[0]).filter(
          (key) => key != 'Samples'
        );

        dispatchMsBurden({ signatureName: nameOptions[0] });
        dispatchMsAssociation({
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
    if (!loadingMsBurden && (fn == 'all' || fn == 'burden')) {
      args.burden = JSON.stringify({
        signatureName: burdenArgs.signatureName,
      });
    }
    if (!loadingMsAssociation && (fn == 'all' || fn == 'association')) {
      args.association = JSON.stringify({
        useCancerType: associationArgs.toggleCancer,
        both: associationArgs.both,
        signatureName1: associationArgs.signatureName1,
        signatureName2: associationArgs.signatureName2,
      });
    }
    if (!loadingMsLandscape && (fn == 'all' || fn == 'landscape')) {
      args.landscape = JSON.stringify({
        variableFile: landscapeArgs.variableFile,
      });
    }
    if (!loadingMsPrevalence && (fn == 'all' || fn == 'prevalence')) {
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
    try {
      const { debugR, output, errors, projectID: pID } = await submitR(
        rFn,
        args,
        id
      );

      if (output) {
        if (!projectID) dispatchExpExposure({ projectID: pID });

        if (fn == 'all') {
          dispatchTMB({
            plotPath: output.tmbPath,
            err: errors.tmbError,
            debugR: debugR,
          });

          dispatchTmbSignatures({
            plotPath: output.signaturePath,
            err: errors.signaturesError,
            debugR: debugR,
          });

          dispatchMsDecomposition({
            plotPath: output.decompositionPath,
            txtPath: output.decompositionData,
            err: errors.decompositionError,
            debugR: debugR,
          });
        }

        if (fn == 'all' || fn == 'burden') {
          dispatchMsBurden({
            plotPath: output.burdenPath,
            err: errors.burdenError,
            debugR: debugR,
          });
        }

        if (fn == 'all' || fn == 'association')
          dispatchMsAssociation({
            plotPath: output.associationPath,
            err: errors.associationError,
            debugR: debugR,
          });

        if (fn == 'all' || fn == 'landscape')
          dispatchMsLandscape({
            plotPath: output.landscapePath,
            err: errors.landscapeError,
            debugR: debugR,
          });

        if (fn == 'all' || fn == 'prevalence')
          dispatchMsPrevalence({
            plotPath: output.prevalencePath,
            err: errors.prevalenceError,
            debugR: debugR,
          });
      } else {
        dispatchError(debugR);
      }
    } catch (err) {
      dispatchError(err);
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

    window.location.hash = '#/exploring/exposure';

    dispatchExpExposure({
      ...initialState.expExposure,
      source: source,
      display: display,
      study: 'PCAWG',
      strategy: 'WGS',
      refSignatureSet: 'COSMIC v3 Signatures (SBS)',
      cancer: 'Lung-AdenoCA',
      studyOptions: studyOptions,
      strategyOptions: strategyOptions,
      refSignatureSetOptions: refSignatureSetOptions,
      cancerOptions: cancerOptions,
      signatureNameOptions: signatureNameOptions,
    });
    dispatchTMB(initialState.tmb);
    dispatchTmbSignatures(initialState.tmbSignatures);
    dispatchMsBurden(initialState.msBurden);
    dispatchMsDecomposition(initialState.msDecomposition);
    dispatchMsAssociation(initialState.msAssociation);
    dispatchMsLandscape(initialState.msLandscape);
    dispatchMsPrevalence(initialState.msPrevalence);
    // populateControls();
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
    dispatchMsLandscape({ variableFile: file.name });
  }

  const tabs = [
    {
      component: <TMB />,
      key: 'tmb',
      name: 'TMB',
      title: 'Tumor Mutational Burden',
    },
    {
      component: <TmbSig />,
      key: 'tmbSig',
      name: 'TMB by Signatures',
      title: 'Tumor Mutational Burden Separated by Signatures',
    },
    {
      component: <MsBurden calculateBurden={calculateBurden} />,
      key: 'msBurden',
      name: 'MS Burden Across Cancer Types',
      title: 'Mutational Signature Burden Across Cancer Types',
    },
    {
      component: <MsDecomposition />,
      key: 'msDecomposition',
      name: 'MS Decomposition',
      title: 'Evaluating the Performance of Mutational Signature Decomposition',
    },
    {
      component: (
        <MsAssociation
          calculateAssociation={calculateAssociation}
          handleSet={handleSet}
          getSignatureNames={getSignatureNames}
        />
      ),
      key: 'msAssociation',
      name: 'MS Association',
      title: 'Mutational Signature Association',
    },
    {
      component: (
        <MsLandscape
          calculateLandscape={calculateLandscape}
          handleVariable={handleVariable}
        />
      ),
      key: 'msLandscape',
      name: 'MS Activity Landscape',
      title: 'Landscape of Mutational Signature Activity',
    },
    {
      component: <MsPrevalence calculatePrevalence={calculatePrevalence} />,
      key: 'msPrevalence',
      name: 'MS Prevalence',
      title: 'Prevalence of Mutational Signature',
    },
  ];

  const examples = [
    {
      title:
        'PCAWG/WGS/COSMIC v3 Signatures (SBS)/ Lung-AdenoCA; MSA SBS5 vs SBS40',
      path: 'exposure1',
    },
    {
      title:
        'PCAWG/WGS/COSMIC v3 Signatures (SBS)/ Skin-Melanoma; MSA SBS7a vs SBS7b',
      path: 'exposure2',
    },
    {
      title:
        'PCAWG/WGS/COSMIC v3 Signatures (SBS)/ Breast-AdenoCA; MSA SBS3 vs SBS5',
      path: 'exposure3',
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
              <Col lg="auto">
                <Group>
                  <Label className="mr-4">Data Source</Label>
                  <Check inline id="radioPublic">
                    <Check.Input
                      disabled={loading || projectID}
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
                      disabled={loading || projectID}
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
                        disabled={loading || projectID}
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
                        disabled={loading || projectID}
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
                        disabled={loading || projectID}
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
                        disabled={loading || projectID}
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
                  <Col>
                    <Group controlId="toggleCancerType">
                      <Check
                        disabled={loading || projectID}
                        type="checkbox"
                        label="Cancer Type Only"
                        value={associationArgs.toggleCancer}
                        checked={associationArgs.toggleCancer}
                        onChange={(e) => {
                          if (!associationArgs.toggleCancer == false)
                            handleSet(refSignatureSet);
                          else {
                            getSignatureNames();
                          }
                          dispatchMsAssociation({
                            toggleCancer: !associationArgs.toggleCancer,
                          });
                        }}
                      />
                    </Group>
                  </Col>
                </Row>
                <Row>
                  <Col lg="6">
                    <Button
                      disabled={loading}
                      className="w-100 mb-3"
                      variant="secondary"
                      onClick={() => handleReset()}
                    >
                      Reset
                    </Button>
                  </Col>
                  <Col lg="6">
                    <Button
                      disabled={loading || projectID}
                      className="w-100"
                      variant="primary"
                      onClick={() => calculateAll()}
                    >
                      Calculate
                    </Button>
                  </Col>
                </Row>
                <hr />
                <strong>Example Queries</strong>
                {examples.map(({ title, external, path }, index) => (
                  <div key={index} className="mb-2">
                    <a href={`#/exploring/exposure/${path}`}>{title}</a>
                    {external && (
                      <span>
                        {'; '}
                        <a href={external.href} target="_blank">
                          {external.name}
                        </a>
                      </span>
                    )}
                  </div>
                ))}
              </div>
            ) : (
              <div>
                <Row>
                  <Col>
                    <Group>
                      <Label>Upload Exposure File</Label>
                      <Form.File
                        disabled={loading || projectID}
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
                        disabled={loading || projectID}
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
                          disabled={loading || projectID}
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
                            disabled={loading || projectID}
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
                            disabled={loading || projectID}
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
                          disabled={loading || projectID}
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
                        disabled={loading || projectID}
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
                  <Col lg="6">
                    <Button
                      disabled={loading}
                      className="w-100 mb-3"
                      variant="secondary"
                      onClick={() => handleReset()}
                    >
                      Reset
                    </Button>
                  </Col>
                  <Col lg="6">
                    <Button
                      disabled={loading || projectID}
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
          <Container
            transition={false}
            className="mt-2"
            defaultActiveKey={display}
            activeKey={display}
            onSelect={(tab) => dispatchExpExposure({ display: tab })}
          >
            <Nav variant="tabs">
              {tabs.map(({ key, name, title }) => (
                <Item key={key}>
                  <Link
                    eventKey={key}
                    as="button"
                    className="outline-none"
                    title={title}
                  >
                    <strong>{name}</strong>
                  </Link>
                </Item>
              ))}
            </Nav>
            <Content
              className={`bg-white tab-pane-bordered rounded-0 d-block`}
              style={{ overflowX: 'auto' }}
            >
              {tabs.map(({ key, component }) => (
                <Pane key={key} eventKey={key} className="border-0">
                  <LoadingOverlay
                    active={loading}
                    content={loadingMsg}
                    showIndicator={loadingMsg}
                  />
                  {component}
                </Pane>
              ))}
            </Content>
          </Container>
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
