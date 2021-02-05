import React, { useState, useEffect, useRef } from 'react';
import { Form, Row, Col, Tab, Button, Nav } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import TMB from './tmb';
import TmbSig from './tmbSignatures';
import MsBurden from './msBurden';
import MsAssociation from './msAssociation';
import MsDecomposition from './msDecomposition';
import MsLandscape from './msLandscape';
import MsPrevalence from './msPrevalence';
import MSIndividual from './msIndividual';
import {
  actions as exploringActions,
  getInitialState,
} from '../../../../services/store/exploring';
import { actions as modalActions } from '../../../../services/store/modal';
import { unique2d } from '../../../../services/utils';
import Select from '../../../controls/select/select';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '../../../controls/sidebar-container/sidebar-container';

const actions = { ...exploringActions, ...modalActions };
const { Group, Label, Check } = Form;
const { Container, Content, Pane } = Tab;
const { Item, Link } = Nav;

export default function Exposure({ match }) {
  const dispatch = useDispatch();
  const exploring = useSelector((state) => state.exploring);

  const mergeExploring = (state) =>
    dispatch(actions.mergeExploring({ exploring: state }));
  const mergeExposure = (state) =>
    dispatch(actions.mergeExploring({ exposure: state }));
  const mergeTMB = (state) => dispatch(actions.mergeExploring({ tmb: state }));
  const mergeTmbSignatures = (state) =>
    dispatch(actions.mergeExploring({ tmbSignatures: state }));
  const mergeMsBurden = (state) =>
    dispatch(actions.mergeExploring({ msBurden: state }));
  const mergeMsAssociation = (state) =>
    dispatch(actions.mergeExploring({ msAssociation: state }));
  const mergeMsDecomposition = (state) =>
    dispatch(actions.mergeExploring({ msDecomposition: state }));
  const mergeMsPrevalence = (state) =>
    dispatch(actions.mergeExploring({ msPrevalence: state }));
  const mergeMsLandscape = (state) =>
    dispatch(actions.mergeExploring({ msLandscape: state }));
  const mergeMsIndividual = (state) =>
    dispatch(actions.mergeExploring({ msIndividual: state }));
  const mergeError = (state) => dispatch(actions.mergeModal({ error: state }));

  const { exampleName } = match.params;

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
    publicSampleOptions,
    userSampleOptions,
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
  } = exploring.exposure;
  const {
    exposureSignature,
    exposureCancer,
    signatureNames,
  } = exploring.exploring;
  const { loading: loadingMsBurden, ...burdenArgs } = exploring.msBurden;
  const {
    loading: loadingMsAssociation,
    ...associationArgs
  } = exploring.msAssociation;
  const {
    loading: loadingMsLandscape,
    ...landscapeArgs
  } = exploring.msLandscape;
  const {
    loading: loadingMsPrevalence,
    ...prevalenceArgs
  } = exploring.msPrevalence;
  const {
    loading: loadingMsIndividual,
    ...individualArgs
  } = exploring.msIndividual;

  const [exposureFileObj, setExposure] = useState(new File([], ''));
  const [matrixFileObj, setMatrix] = useState(new File([], ''));
  const [signatureFileObj, setSignature] = useState(new File([], ''));
  const [variableFileObj, setVariable] = useState(new File([], ''));

  // const [exposureValidity, setExposureValidity] = useState(false);
  // const [matrixValidity, setMatrixValidity] = useState(false);
  // const [signatureValidity, setSignatureValidity] = useState(false);

  // load example if available
  useEffect(() => {
    if (exampleName) loadExample(exampleName);
  }, [exampleName]);

  // set selected tab on component render
  useEffect(() => {
    mergeExploring({ displayTab: 'exposure' });
  }, []);

  function usePrevious(value) {
    const ref = useRef();
    useEffect(() => {
      ref.current = value;
    }, [value]);

    return ref.current;
  }

  // get new signature name options filtered by cancer type on first mount
  const prevSigNameOptions = usePrevious(signatureNameOptions || []);
  useEffect(() => {
    if (
      source == 'public' &&
      associationArgs.toggleCancer &&
      study &&
      strategy &&
      refSignatureSet &&
      prevSigNameOptions
    )
      getSignatureNames();
  }, [study, strategy, refSignatureSet, cancer, associationArgs.toggleCancer]);

  // get sample name options
  const prevSampleOptions = usePrevious(publicSampleOptions);
  useEffect(() => {
    if (
      source == 'public' &&
      study &&
      strategy &&
      refSignatureSet &&
      cancer &&
      prevSampleOptions
    )
      getSampleNames();
  }, [study, strategy, refSignatureSet, cancer]);

  async function loadExample(id) {
    mergeExposure({
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

      mergeExposure({ ...state.expExposure, projectID: projectID });
      // rehydrate state if available
      if (state.tmb) mergeTMB(state.tmb);
      if (state.msBurden) mergeMsBurden(state.msBurden);
      if (state.msAssociation) mergeMsAssociation(state.msAssociation);
      if (state.msDecomposition) mergeMsDecomposition(state.msDecomposition);
      if (state.msLandscape) mergeMsLandscape(state.msLandscape);
      if (state.msPrevalence) mergeMsPrevalence(state.msPrevalence);
      if (state.tmbSignatures) mergeTmbSignatures(state.tmbSignatures);
    } catch (error) {
      mergeError({ visible: true, message: error.message });
    }
    mergeExposure({
      loading: false,
    });
  }

  // get signature name options filtered by cancer type
  async function getSignatureNames() {
    mergeExposure({
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
        mergeExposure({
          signatureNameOptions: output.data,
        });
      else mergeError(stdout);
    } catch (err) {
      mergeError({ visible: true, message: err.message });
    }
    mergeExposure({ loading: false, loadingMsg: null });
  }

  // get sample name options filtered by cancer type
  async function getSampleNames() {
    mergeExposure({
      loading: true,
      loadingMsg: 'Filtering Sample Names',
    });
    try {
      const { stdout, output } = await (
        await fetch(`api/getSampleNames`, {
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

      if (output.data.length) {
        mergeExposure({
          publicSampleOptions: output.data,
        });
        mergeMsIndividual({ sample: output.data[0] });
      } else mergeError(stdout);
    } catch (err) {
      mergeError({ visible: true, message: err.message });
    }
    mergeExposure({ loading: false, loadingMsg: null });
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
        mergeError('R submit failed');
      }
    } catch (err) {
      mergeError({ visible: true, message: err.message });
    }
  }

  async function calculateBurden() {
    try {
      if (source == 'user' && !projectID) {
        mergeError('Missing Required Files');
      } else {
        mergeMsBurden({
          loading: true,
          err: false,
          plotPath: '',
        });

        if (source == 'user') {
          if (!projectID) {
            const id = handleUpload();
            await handleCalculate('burden', id);
          } else {
            await handleCalculate('burden');
          }
        } else {
          await handleCalculate('burden');
        }

        mergeMsBurden({ loading: false });
      }
    } catch (error) {
      mergeError({ visible: true, message: error.message });
    }
  }

  async function calculateAssociation() {
    try {
      if (source == 'user' && !projectID) {
        mergeError('Missing Required Files');
      } else {
        mergeMsAssociation({
          loading: true,
          err: false,
          plotPath: '',
        });

        if (source == 'user') {
          if (!projectID) {
            const id = await handleUpload();
            await handleCalculate('association', id);
          } else {
            await handleCalculate('association');
          }
        } else {
          await handleCalculate('association');
        }

        mergeMsAssociation({ loading: false });
      }
    } catch (error) {
      mergeError({ visible: true, message: error.message });
    }
  }

  async function calculateLandscape() {
    try {
      if (source == 'user' && !projectID) {
        mergeError('Missing Required Files');
      } else {
        mergeMsLandscape({
          loading: true,
          err: false,
          plotPath: '',
        });

        if (source == 'user') {
          if (!projectID) {
            const id = await handleUpload();
            await handleCalculate('landscape', id);
          } else {
            await handleCalculate('landscape');
          }
        } else if (variableFileObj.size) {
          const id = await uploadVariable();
          await handleCalculate('landscape', id);
        } else {
          await handleCalculate('landscape');
        }

        mergeMsLandscape({ loading: false });
      }
    } catch (error) {
      mergeError({ visible: true, message: error.message });
    }
  }

  async function calculatePrevalence() {
    try {
      if (source == 'user' && !projectID) {
        mergeError('Missing Required Files');
      } else {
        mergeMsPrevalence({
          loading: true,
          err: false,
          plotPath: '',
        });

        if (source == 'user') {
          if (!projectID) {
            const id = await handleUpload();
            await handleCalculate('prevalence', id);
          } else {
            await handleCalculate('prevalence');
          }
        } else {
          await handleCalculate('prevalence');
        }

        mergeMsPrevalence({ loading: false });
      }
    } catch (error) {
      mergeError({ visible: true, message: error.message });
    }
  }
  async function calculateIndividual() {
    try {
      if (source == 'user' && !projectID) {
        mergeError('Missing Required Files');
      } else {
        mergeMsIndividual({
          loading: true,
          err: false,
          plotPath: '',
        });

        if (source == 'user') {
          if (!projectID) {
            const id = await handleUpload();
            await handleCalculate('individual', id);
          } else {
            await handleCalculate('individual');
          }
        } else {
          await handleCalculate('individual');
        }

        mergeMsIndividual({ loading: false });
      }
    } catch (error) {
      mergeError({ visible: true, message: error.message });
    }
  }

  async function calculateAll() {
    try {
      mergeExposure({ loading: true });

      mergeTMB({
        plotPath: '',
        err: '',
      });
      mergeTmbSignatures({
        plotPath: '',
        err: '',
      });
      mergeMsDecomposition({
        plotPath: '',
        txtPath: '',
        err: '',
      });
      mergeMsBurden({
        plotPath: '',
        err: '',
      });
      mergeMsAssociation({
        plotPath: '',
        err: '',
      });
      mergeMsLandscape({
        plotPath: '',
        err: '',
      });

      mergeMsPrevalence({
        plotPath: '',
        err: '',
      });

      if (source == 'user') {
        const { projectID, exposureData } = await handleUpload();

        // get signature name options, ignore sample key
        const nameOptions = exposureData.columns.filter(
          (key) => key != 'Samples'
        );

        const sampleOptions = unique2d(
          'Samples',
          exposureData.columns,
          exposureData.data
        );

        mergeMsBurden({ signatureName: nameOptions[0] });
        mergeMsAssociation({
          signatureName1: nameOptions[0],
          signatureName2: nameOptions[1],
        });
        mergeMsIndividual({ sample: sampleOptions[0] });
        mergeExposure({
          projectID: projectID,
          userNameOptions: nameOptions,
          userSampleOptions: sampleOptions,
        });

        await handleCalculate('all', projectID);
      } else if (variableFileObj.size) {
        const id = await uploadVariable();
        await handleCalculate('all', id);
      } else {
        await handleCalculate('all');
      }
    } catch (err) {
      mergeError({ visible: true, message: err.message });
    } finally {
      mergeExposure({ loading: false });
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
    if (!loadingMsIndividual && (fn == 'all' || fn == 'individual')) {
      args.individual = JSON.stringify({
        sample: individualArgs.sample,
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
        if (!projectID) mergeExposure({ projectID: pID });

        if (fn == 'all') {
          mergeTMB({
            plotPath: output.tmbPath,
            err: errors.tmbError,
            debugR: debugR,
          });

          mergeTmbSignatures({
            plotPath: output.signaturePath,
            err: errors.signaturesError,
            debugR: debugR,
          });

          mergeMsDecomposition({
            plotPath: output.decompositionPath,
            txtPath: output.decompositionData,
            err: errors.decompositionError,
            debugR: debugR,
          });
        }

        if (fn == 'all' || fn == 'burden') {
          mergeMsBurden({
            plotPath: output.burdenPath,
            err: errors.burdenError,
            debugR: debugR,
          });
        }

        if (fn == 'all' || fn == 'association')
          mergeMsAssociation({
            plotPath: output.associationPath,
            err: errors.associationError,
            debugR: debugR,
          });

        if (fn == 'all' || fn == 'landscape')
          mergeMsLandscape({
            plotPath: output.landscapePath,
            err: errors.landscapeError,
            debugR: debugR,
          });

        if (fn == 'all' || fn == 'prevalence')
          mergeMsPrevalence({
            plotPath: output.prevalencePath,
            err: errors.prevalenceError,
            debugR: debugR,
          });

        if (fn == 'all' || fn == 'individual')
          mergeMsIndividual({
            plotPath: output.individualPath,
            err: errors.individualError,
            debugR: debugR,
          });
      } else {
        mergeError(debugR);
      }
    } catch (err) {
      mergeError({ visible: true, message: err.message });
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

    mergeExposure({
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

    mergeExposure({
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

    mergeExposure({
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
            mergeError(message);
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
          mergeError({ visible: true, message: err.message });
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

    mergeExposure({
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
      publicSampleOptions: publicSampleOptions,
    });
    mergeTMB(initialState.tmb);
    mergeTmbSignatures(initialState.tmbSignatures);
    mergeMsBurden(initialState.msBurden);
    mergeMsDecomposition(initialState.msDecomposition);
    mergeMsAssociation(initialState.msAssociation);
    mergeMsLandscape(initialState.msLandscape);
    mergeMsPrevalence(initialState.msPrevalence);
    mergeMsIndividual(initialState.msIndividual);
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
          mergeError(message);
          reject(error);
        } else {
          const { projectID } = await response.json();
          resolve(projectID);
        }
      } catch (err) {
        mergeError({ visible: true, message: err.message });
        reject(err);
      }
    });
  }

  function handleVariable(file) {
    setVariable(file);
    mergeMsLandscape({ variableFile: file.name });
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
      name: 'TMB Signatures',
      title: 'Tumor Mutational Burden Separated by Signatures',
    },
    {
      component: <MsBurden calculateBurden={calculateBurden} />,
      key: 'msBurden',
      name: 'MS Burden',
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
      name: 'MS Landscape',
      title: 'Landscape of Mutational Signature Activity',
    },
    {
      component: <MsPrevalence calculatePrevalence={calculatePrevalence} />,
      key: 'msPrevalence',
      name: 'MS Prevalence',
      title: 'Prevalence of Mutational Signature',
    },
    {
      component: <MSIndividual calculateIndividual={calculateIndividual} />,
      key: 'msIndividaul',
      name: 'MS Individual',
      title: 'Mutational Signature in Individual Sample',
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
        onCollapsed={(e) => mergeExposure({ openSidebar: !e })}
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
                      onChange={(e) => mergeExposure({ source: 'public' })}
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
                      onChange={(e) => mergeExposure({ source: 'user' })}
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
                          mergeExposure({
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
                          mergeMsAssociation({
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
                          mergeExposure({
                            exposureFile: e.target.files[0].name,
                          });
                        }}
                        custom
                      />
                      {/* {exposureValidity && (
                        <span className="text-danger">
                          Exposure File Required
                        </span>
                      )} */}
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
                          mergeExposure({
                            matrixFile: e.target.files[0].name,
                          });
                        }}
                        custom
                      />
                      {/* {matrixValidity && (
                        <span className="text-danger">
                          Matrix File Required
                        </span>
                      )} */}
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
                            mergeExposure({
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
                            mergeExposure({
                              signatureFile: e.target.files[0].name,
                            });
                          }}
                          custom
                        />
                        {/* {signatureValidity && (
                          <span className="text-danger">
                            Signature File Required
                          </span>
                        )} */}
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
                        onChange={(genome) => mergeExposure({ genome: genome })}
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
            onSelect={(tab) => mergeExposure({ display: tab })}
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
