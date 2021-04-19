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
  actions as explorationActions,
  getInitialState,
} from '../../../../services/store/exploration';
import { actions as modalActions } from '../../../../services/store/modal';
import { unique2d } from '../../../../services/utils';
import Select from '../../../controls/select/select';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '../../../controls/sidebar-container/sidebar-container';

const actions = { ...explorationActions, ...modalActions };
const { Group, Label, Check } = Form;
const { Container, Content, Pane } = Tab;
const { Item, Link } = Nav;

export default function Exposure({ match }) {
  const dispatch = useDispatch();
  const exploration = useSelector((state) => state.exploration);

  const mergeState = async (state) =>
    await dispatch(actions.mergeExploration({ ...state }));
  const mergeExposure = (state) =>
    dispatch(actions.mergeExploration({ exposure: state }));
  const mergeTMB = (state) =>
    dispatch(actions.mergeExploration({ tmb: state }));
  const mergeTmbSignatures = (state) =>
    dispatch(actions.mergeExploration({ tmbSignatures: state }));
  const mergeMsBurden = (state) =>
    dispatch(actions.mergeExploration({ msBurden: state }));
  const mergeMsAssociation = (state) =>
    dispatch(actions.mergeExploration({ msAssociation: state }));
  const mergeMsDecomposition = (state) =>
    dispatch(actions.mergeExploration({ msDecomposition: state }));
  const mergeMsPrevalence = (state) =>
    dispatch(actions.mergeExploration({ msPrevalence: state }));
  const mergeMsLandscape = (state) =>
    dispatch(actions.mergeExploration({ msLandscape: state }));
  const mergeMsIndividual = (state) =>
    dispatch(actions.mergeExploration({ msIndividual: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const { exampleName } = match.params;

  const {
    study,
    studyOptions,
    strategy,
    strategyOptions,
    cancer,
    cancerOptions,
    refSignatureSet,
    refSignatureSetOptions,
    useCancerType,
    signatureNameOptions,
    publicSampleOptions,
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
    submitted,
    submitAll,
    gettingSignatureNames,
    gettingSampleNames,
  } = exploration.exposure;
  const {
    exposureSignature,
    exposureCancer,
    signatureNames,
  } = exploration.exploration;

  const { loading: tmbLoading } = exploration.tmb;
  const { loading: tmbSigLoading } = exploration.tmbSignatures;
  const { loading: decompositionLoading } = exploration.msDecomposition;
  const { loading: loadingMsBurden, ...burdenArgs } = exploration.msBurden;
  const {
    loading: loadingMsAssociation,
    ...associationArgs
  } = exploration.msAssociation;
  const {
    loading: loadingMsLandscape,
    ...landscapeArgs
  } = exploration.msLandscape;
  const {
    loading: loadingMsPrevalence,
    ...prevalenceArgs
  } = exploration.msPrevalence;
  const {
    loading: loadingMsIndividual,
    ...individualArgs
  } = exploration.msIndividual;

  const [exposureFileObj, setExposure] = useState(new File([], ''));
  const [matrixFileObj, setMatrix] = useState(new File([], ''));
  const [signatureFileObj, setSignature] = useState(new File([], ''));
  const [variableFileObj, setVariable] = useState(new File([], ''));

  const [exposureValidity, setExposureValidity] = useState(false);
  const [matrixValidity, setMatrixValidity] = useState(false);
  const [signatureValidity, setSignatureValidity] = useState(false);

  const [checkValid, setCheckValid] = useState(false);
  const [expand, setExpand] = useState(false);

  // load example if available
  useEffect(() => {
    if (exampleName) loadExample(exampleName);
  }, [exampleName]);

  // set selected tab on component render
  // get new signature name options filtered by cancer type on first render
  useEffect(() => {
    mergeState({ exploration: { displayTab: 'exposure' } });
    if (source == 'public') {
      if (!gettingSignatureNames) getSignatureNames();
      if (!gettingSampleNames) getSampleNames();
    }
  }, []);

  // lazy load plots after loading signature names and sample names filtered by cancer type
  useEffect(() => {
    if (submitAll) {
      if (!gettingSignatureNames) {
        if (
          !loadingMsBurden &&
          !burdenArgs.plotPath &&
          burdenArgs.signatureName
        )
          calculateBurden();
        if (
          !loadingMsAssociation &&
          !associationArgs.plotPath &&
          associationArgs.signatureName1
        )
          calculateAssociation();
      }
      if (
        !gettingSampleNames &&
        !loadingMsIndividual &&
        !individualArgs.plotPath &&
        individualArgs.sample
      )
        calculateIndividual();
    }
  }, [gettingSignatureNames, gettingSampleNames]);

  function usePrevious(value) {
    const ref = useRef();
    useEffect(() => {
      ref.current = value;
    }, [value]);

    return ref.current;
  }

  // get signature name options
  const prevSigNameOptions = usePrevious(signatureNameOptions || []);
  useEffect(() => {
    if (
      source == 'public' &&
      useCancerType &&
      study &&
      strategy &&
      refSignatureSet &&
      prevSigNameOptions
    )
      getSignatureNames();
  }, [study, strategy, refSignatureSet, cancer, useCancerType]);

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
    mergeState({
      exposure: {
        loading: {
          active: true,
          content: 'Loading Example',
          showIndicator: true,
        },
      },
    });
    try {
      const { projectID, state } = await (
        await fetch(`api/getExposureExample/${id}`)
      ).json();

      const { expExposure, ...rest } = state;

      mergeState({
        exposure: { expExposure, projectID: projectID },
        ...rest,
      });
    } catch (error) {
      mergeError(error.message);
    }
    mergeState({
      exposure: {
        loading: false,
      },
    });
  }

  function validateFiles() {
    setCheckValid(true);
    exposureFileObj.size
      ? setExposureValidity(true)
      : setExposureValidity(false);
    matrixFileObj.size ? setMatrixValidity(true) : setMatrixValidity(false);
    signatureFileObj.size
      ? setSignatureValidity(true)
      : setSignatureValidity(false);

    return exposureValidity && matrixValidity && signatureValidity;
  }

  // get signature name options filtered by cancer type
  async function getSignatureNames() {
    mergeState({ exposure: { gettingSignatureNames: true } });
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
        mergeState({
          exposure: {
            signatureNameOptions: output.data,
          },
        });
      else console.log('No Signature Names Found');
    } catch (err) {
      mergeError(err.message);
    }
    mergeState({ exposure: { gettingSignatureNames: false } });
  }

  // get sample name options filtered by cancer type
  async function getSampleNames() {
    mergeState({ exposure: { gettingSampleNames: true } });
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
      } else console.log('No Sample Names Found');
    } catch (err) {
      mergeError(err.message);
    }
    mergeState({ exposure: { gettingSampleNames: false } });
  }

  async function submitR(fn, args, id = projectID) {
    try {
      const response = await fetch(`api/explorationR`, {
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
      mergeError(err.message);
    }
  }

  async function calculateTMB() {
    try {
      if (source == 'user' && !projectID) {
        mergeError('Missing Required Files');
      } else {
        mergeTMB({
          loading: true,
          err: false,
          plotPath: '',
        });

        if (source == 'user') {
          if (!projectID) {
            const id = handleUpload();
            await handleCalculate('tmb', id);
          } else {
            await handleCalculate('tmb');
          }
        } else {
          await handleCalculate('tmb');
        }

        mergeTMB({ loading: false });
      }
    } catch (error) {
      mergeError(error.message);
    }
  }

  async function calculateTmbSig() {
    try {
      if (source == 'user' && !projectID) {
        mergeError('Missing Required Files');
      } else {
        mergeTmbSignatures({
          loading: true,
          err: false,
          plotPath: '',
        });

        if (source == 'user') {
          if (!projectID) {
            const id = handleUpload();
            await handleCalculate('tmbSig', id);
          } else {
            await handleCalculate('tmbSig');
          }
        } else {
          await handleCalculate('tmbSig');
        }

        mergeTmbSignatures({ loading: false });
      }
    } catch (error) {
      mergeError(error.message);
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
      mergeError(error.message);
    }
  }

  async function calculateDecomposition() {
    try {
      if (source == 'user' && !projectID) {
        mergeError('Missing Required Files');
      } else {
        mergeMsDecomposition({
          loading: true,
          err: false,
          plotPath: '',
        });

        if (source == 'user') {
          if (!projectID) {
            const id = handleUpload();
            await handleCalculate('decomposition', id);
          } else {
            await handleCalculate('decomposition');
          }
        } else {
          await handleCalculate('decomposition');
        }

        mergeMsDecomposition({ loading: false });
      }
    } catch (error) {
      mergeError(error.message);
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
      mergeError(error.message);
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
      mergeError(error.message);
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
      mergeError(error.message);
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
      mergeError(error.message);
    }
  }

  function handleLoading(active) {
    mergeTMB({ loading: active });
    mergeTmbSignatures({ loading: active });
    mergeMsDecomposition({ loading: active });
    if (!gettingSignatureNames) mergeMsBurden({ loading: active });
    if (!gettingSignatureNames) mergeMsAssociation({ loading: active });
    mergeMsLandscape({ loading: active });
    mergeMsPrevalence({ loading: active });
    if (!gettingSampleNames) mergeMsIndividual({ loading: active });
  }

  async function calculateAll() {
    try {
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
        const params = {
          exposure: {
            projectID: projectID,
            userNameOptions: nameOptions,
            userSampleOptions: sampleOptions,
          },
          msIndividual: { sample: sampleOptions[0] },
          msAssociation: {
            signatureName1: nameOptions[0],
            signatureName2: nameOptions[1],
          },
          msBurden: { signatureName: nameOptions[0] },
        };
        mergeState(params);

        await handleCalculate('all', projectID, params);
      } else if (variableFileObj.size) {
        const id = await uploadVariable();
        await handleCalculate('all', id);
      } else {
        await handleCalculate('all');
      }
    } catch (err) {
      mergeError(err.message);
    }
  }

  async function handleCalculate(fn = 'all', id = projectID, params = {}) {
    let rFn = 'exposurePublic';
    let args = {
      fn: fn,
      common: JSON.stringify({
        study: study,
        strategy: strategy,
        refSignatureSet: refSignatureSet,
        cancerType: cancer,
        genome: genome,
        useCancerType: useCancerType,
      }),
    };
    if (
      !loadingMsBurden &&
      !gettingSignatureNames &&
      (fn == 'all' || fn == 'burden')
    ) {
      args.burden = JSON.stringify({
        signatureName: burdenArgs.signatureName,
        ...params.msBurden,
      });
    }
    if (
      !loadingMsAssociation &&
      !gettingSignatureNames &&
      (fn == 'all' || fn == 'association')
    ) {
      args.association = JSON.stringify({
        // useCancerType: associationArgs.toggleCancer,
        both: associationArgs.both,
        signatureName1: associationArgs.signatureName1,
        signatureName2: associationArgs.signatureName2,
        ...params.msAssociation,
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
    if (
      !loadingMsIndividual &&
      !gettingSampleNames &&
      (fn == 'all' || fn == 'individual')
    ) {
      args.individual = JSON.stringify({
        sample: individualArgs.sample,
        ...params.msIndividual,
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
      if (fn == 'all') {
        handleLoading(true);
        mergeExposure({ submitAll: true });
      }
      mergeExposure({ submitted: fn });
      const { debugR, output, errors, projectID: pID } = await submitR(
        rFn,
        args,
        id
      );

      if (output) {
        mergeExposure({ projectID: pID });

        if (output.tmbPath) {
          mergeTMB({
            plotPath: output.tmbPath,
            err: errors.tmbError,
            debugR: debugR,
          });
        }

        if (output.signaturePath) {
          mergeTmbSignatures({
            plotPath: output.signaturePath,
            err: errors.signaturesError,
            debugR: debugR,
          });
        }

        if (output.decompositionPath) {
          mergeMsDecomposition({
            plotPath: output.decompositionPath,
            txtPath: output.decompositionData,
            err: errors.decompositionError,
            debugR: debugR,
          });
        }

        if (output.burdenPath) {
          mergeMsBurden({
            plotPath: output.burdenPath,
            err: errors.burdenError,
            debugR: debugR,
          });
        }

        if (output.associationPath)
          mergeMsAssociation({
            plotPath: output.associationPath,
            err: errors.associationError,
            debugR: debugR,
          });

        if (output.landscapePath)
          mergeMsLandscape({
            plotPath: output.landscapePath,
            err: errors.landscapeError,
            debugR: debugR,
          });

        if (output.prevalencePath)
          mergeMsPrevalence({
            plotPath: output.prevalencePath,
            err: errors.prevalenceError,
            debugR: debugR,
          });

        if (output.individualPath)
          mergeMsIndividual({
            plotPath: output.individualPath,
            err: errors.individualError,
            debugR: debugR,
          });
      } else {
        mergeError(debugR);
      }
    } catch (err) {
      mergeError(err.message);
    } finally {
      if (fn == 'all') handleLoading(false);
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
        exposureValidity &&
        matrixValidity &&
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
          mergeError(err.message);
          reject(err);
        }
      } else {
        reject('Missing required files');
      }
    });
  }

  function handleReset() {
    const initialState = getInitialState();
    setCheckValid(false);
    window.location.hash = '#/exploration/exposure';

    mergeExposure({
      ...initialState.exposure,
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
        mergeError(err.message);
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
      component: <TMB calculateTMB={calculateTMB} />,
      key: 'tmb',
      name: 'TMB',
      title: 'Tumor Mutational Burden',
    },
    {
      component: <TmbSig calculateTmbSig={calculateTmbSig} />,
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
      component: (
        <MsDecomposition calculateDecomposition={calculateDecomposition} />
      ),
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

  const queries = [
    {
      study: 'PCAWG',
      examples: [
        {
          name: 'Lung',
          title:
            'PCAWG/WGS/COSMIC v3 Signatures (SBS)/ Lung-AdenoCA; MSA SBS5 vs SBS40',
          path: 'exposure1',
        },
        {
          name: 'Skin',
          title:
            'PCAWG/WGS/COSMIC v3 Signatures (SBS)/ Skin-Melanoma; MSA SBS7a vs SBS7b',
          path: 'exposure2',
        },
        {
          name: 'Breast',
          title:
            'PCAWG/WGS/COSMIC v3 Signatures (SBS)/ Breast-AdenoCA; MSA SBS3 vs SBS5',
          path: 'exposure3',
        },
      ],
    },
    {
      study: 'TBA',
      examples: [],
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
            <strong>Example Queries</strong>
            {queries.map(({ study, examples }, index) => {
              return index == 0 ? (
                <Row key={`${study}-${index}`} className="mb-2">
                  <Col md="3">{study}:</Col>
                  {examples.map(({ name, title, path }) => (
                    <Col md="3" key={name + index}>
                      <span className="mb-2" disabled={loading || submitted}>
                        <a
                          href={`#/exploration/exposure/${path}`}
                          title={title}
                        >
                          {name}
                        </a>
                      </span>
                    </Col>
                  ))}
                </Row>
              ) : (
                <>
                  {expand ? (
                    <>
                      <Row className="mb-2">
                        <Col md="3">{study}:</Col>
                        {examples.map(({ name, title, path }) => (
                          <Col md="3" key={`${study}-${name}`}>
                            <span className="mb-2">
                              <a
                                href={`#/exploration/exposure/${path}`}
                                title={title}
                              >
                                {name}
                              </a>
                            </span>
                          </Col>
                        ))}
                      </Row>
                      <Row>
                        <Col md="6">
                          <Button
                            onClick={() => setExpand(false)}
                            variant="link"
                            className="p-0"
                            style={{ textDecoration: 'none' }}
                          >
                            Show Less
                          </Button>
                        </Col>
                      </Row>
                    </>
                  ) : (
                    <Row>
                      <Col md="6">
                        <Button
                          onClick={() => setExpand(true)}
                          variant="link"
                          className="p-0"
                          style={{ textDecoration: 'none' }}
                        >
                          Show More
                        </Button>
                      </Col>
                    </Row>
                  )}
                </>
              );
            })}
            <hr className="mb-2" />
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
                        disabled={
                          loading ||
                          submitted ||
                          gettingSignatureNames ||
                          gettingSampleNames
                        }
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
                        disabled={
                          loading ||
                          submitted ||
                          gettingSignatureNames ||
                          gettingSampleNames
                        }
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
                        disabled={
                          loading ||
                          submitted ||
                          gettingSignatureNames ||
                          gettingSampleNames
                        }
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
                        disabled={
                          loading ||
                          submitted ||
                          gettingSignatureNames ||
                          gettingSampleNames
                        }
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
                        disabled={
                          loading ||
                          submitted ||
                          gettingSampleNames ||
                          gettingSignatureNames
                        }
                        type="checkbox"
                        label="Cancer Type Only"
                        value={useCancerType}
                        checked={useCancerType}
                        onChange={(e) => {
                          if (!useCancerType == false)
                            handleSet(refSignatureSet);
                          else {
                            getSignatureNames();
                          }
                          mergeExposure({
                            useCancerType: !useCancerType,
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
                      disabled={loading || submitted}
                      className="w-100"
                      variant="primary"
                      onClick={() => calculateAll()}
                    >
                      Calculate All
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
                        disabled={loading || submitted}
                        id="uploadExposure"
                        label={exposureFileObj.name || 'Exposure File'}
                        accept=".txt"
                        isInvalid={checkValid ? !exposureValidity : false}
                        feedback="Upload an exposure file"
                        onChange={(e) => {
                          setExposure(e.target.files[0]);
                          mergeExposure({
                            exposureFile: e.target.files[0].name,
                          });
                        }}
                        custom
                      />
                    </Group>
                  </Col>
                </Row>
                <Row>
                  <Col>
                    <Group>
                      <Label>Upload Matrix File</Label>
                      <Form.File
                        disabled={loading || submitted}
                        id="uploadMatrix"
                        label={matrixFileObj.name || 'Matrix File'}
                        accept=".txt"
                        isInvalid={checkValid ? !matrixValidity : false}
                        feedback="Upload a matrix file"
                        onChange={(e) => {
                          setMatrix(e.target.files[0]);
                          mergeExposure({
                            matrixFile: e.target.files[0].name,
                          });
                        }}
                        custom
                      />
                    </Group>
                  </Col>
                </Row>
                <Row>
                  <Col>
                    <Group controlId="toggleSignatureSource" className="d-flex">
                      <Label className="mr-4">Use Public Signature Data</Label>
                      <Check inline id="toggleSignatureSource">
                        <Check.Input
                          disabled={loading || submitted}
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
                            disabled={loading || submitted}
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
                            disabled={loading || submitted}
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
                          disabled={loading || submitted}
                          id="uploadSignature"
                          label={signatureFileObj.name || 'Signature File'}
                          accept=".txt"
                          isInvalid={checkValid ? !signatureValidity : false}
                          feedback="Upload a signature file"
                          onChange={(e) => {
                            setSignature(e.target.files[0]);
                            mergeExposure({
                              signatureFile: e.target.files[0].name,
                            });
                          }}
                          custom
                        />
                      </Group>
                    </Col>
                  </Row>
                )}
                <Row>
                  <Col>
                    <Group>
                      <Select
                        disabled={loading || submitted}
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
                      disabled={loading}
                      className="w-100"
                      variant="primary"
                      onClick={() => {
                        if (validateFiles()) calculateAll();
                      }}
                    >
                      Calculate All
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
                <Pane
                  key={key}
                  eventKey={key}
                  className="border-0"
                  style={{ minHeight: '7rem' }}
                >
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
