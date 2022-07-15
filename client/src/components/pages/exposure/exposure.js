import React, { useState, useEffect } from 'react';
import { Form, Row, Col, Button, Nav } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { saveAs } from 'file-saver';
import PublicForm from './publicForm/publicForm';
// import PublicForm from './publicForm';
import UserForm from './userForm';
import Instructions from './instructions';
import TMB from './tmb';
import TMB2 from './tmb/tmb.js';
import TmbSig from './tmbSignatures';
import TmbSig2 from './tmbSignature/tmbSignature.js';
import MsBurden from './msBurden';
import MsBurden2 from './msBurden/msBurden.js';
import MsAssociation from './msAssociation';
import MsDecomposition from './msDecomposition';
import MsLandscape from './msLandscape';
import MsPrevalence from './msPrevalence';
import MSIndividual from './msIndividual';
import Download from './download';
import { actions as exposureActions } from '../../../services/store/exposure';
import { actions as modalActions } from '../../../services/store/modal';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '../../controls/sidebar-container/sidebar-container';

const actions = { ...exposureActions, ...modalActions };
const { Group, Label, Check } = Form;

export default function Exposure({ match }) {
  const dispatch = useDispatch();

  const mergeState = (state) =>
    dispatch(actions.mergeExposure({ main: state }));
  const mergeTMB = (state) => dispatch(actions.mergeExposure({ tmb: state }));
  const mergeTmbSignatures = (state) =>
    dispatch(actions.mergeExposure({ tmbSignatures: state }));
  const mergeMsBurden = (state) =>
    dispatch(actions.mergeExposure({ msBurden: state }));
  const mergeMsAssociation = (state) =>
    dispatch(actions.mergeExposure({ msAssociation: state }));
  const mergeMsDecomposition = (state) =>
    dispatch(actions.mergeExposure({ msDecomposition: state }));
  const mergeMsPrevalence = (state) =>
    dispatch(actions.mergeExposure({ msPrevalence: state }));
  const mergeMsLandscape = (state) =>
    dispatch(actions.mergeExposure({ msLandscape: state }));
  const mergeMsIndividual = (state) =>
    dispatch(actions.mergeExposure({ msIndividual: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));
  const resetExposure = (_) => dispatch(actions.resetExposure());

  const exposureStore = useSelector((state) => state.exposure);

  const { exampleName } = match.params;

  const [variableFileObj, setVariable] = useState(new File([], ''));

  const {
    displayTab,
    exposureSignature,
    exposureCancer,
    study,
    studyOptions,
    strategy,
    strategyOptions,
    cancer,
    cancerOptions,
    rsSet,
    rsSetOptions,
    useCancerType,
    genome,
    exposureFile,
    matrixFile,
    signatureFile,
    source,
    loading,
    projectID,
    openSidebar,
    submitted,
  } = exposureStore.main;

  const { loading: loadingMsBurden, ...burdenArgs } = exposureStore.msBurden;
  const { loading: loadingMsAssociation, ...associationArgs } =
    exposureStore.msAssociation;
  const { loading: loadingMsLandscape, ...landscapeArgs } =
    exposureStore.msLandscape;
  const { loading: loadingMsPrevalence, ...prevalenceArgs } =
    exposureStore.msPrevalence;
  const { loading: loadingMsIndividual, ...individualArgs } =
    exposureStore.msIndividual;

  // load example if available
  useEffect(() => {
    if (exampleName) loadExample(exampleName);
  }, [exampleName]);

  async function loadExample(id) {
    mergeState({
      loading: {
        active: true,
        // content: 'Loading Example',
        // showIndicator: true,
      },
    });
    try {
      const { state } = await (
        await fetch(`web/getExposureExample/${id}`)
      ).json();

      dispatch(actions.mergeExposure(state));
    } catch (error) {
      mergeError('Example does not exist');
    }
    mergeState({
      loading: false,
      submitted: true,
      openSidebar: false,
    });
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

        await handleCalculate('burden');

        mergeMsBurden({ loading: false });
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

        await handleCalculate('association');

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
        const id = variableFileObj ? await uploadVariable() : projectID;
        console.log(variableFileObj);
        await handleCalculate('landscape', id);
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

        await handleCalculate('prevalence');

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

        await handleCalculate('individual');

        mergeMsIndividual({ loading: false });
      }
    } catch (error) {
      mergeError(error.message);
    }
  }

  async function submitR(fn, args, id = projectID) {
    try {
      const response = await fetch(`web/explorationWrapper`, {
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

  async function handleCalculate(fn = 'all', id = projectID, params = {}) {
    mergeState({ loading: true, submitted: true });

    let rFn = 'exposurePublic';
    let args = {
      fn: fn,
      common: JSON.stringify({
        study: study,
        strategy: strategy,
        rsSet: rsSet,
        cancerType: cancer,
        genome: genome,
        useCancerType: useCancerType,
      }),
    };
    if (fn == 'all' || fn == 'burden') {
      args.burden = JSON.stringify({
        signatureName: burdenArgs.signatureName,
      });
    }
    if (fn == 'all' || fn == 'association') {
      args.association = JSON.stringify({
        both: associationArgs.both,
        signatureName1: associationArgs.signatureName1,
        signatureName2: associationArgs.signatureName2,
      });
    }
    if (fn == 'all' || fn == 'landscape') {
      args.landscape = JSON.stringify({
        variableFile: landscapeArgs.variableFile,
      });
    }
    if (fn == 'all' || fn == 'prevalence') {
      args.prevalence = JSON.stringify({
        mutation: parseFloat(prevalenceArgs.mutation) || 100,
      });
    }
    if (fn == 'all' || fn == 'individual') {
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
      const { stdout, output, projectID: pID } = await submitR(rFn, args, id);

      if (Object.keys(output).length) {
        mergeState({ projectID: pID });

        mergeTMB({
          plotPath: output.tmbPath,
          err: output.tmbError,
        });

        mergeTmbSignatures({
          plotPath: output.signaturePath,
          err: output.signaturesError,
        });

        mergeMsDecomposition({
          plotPath: output.decompositionPath,
          txtPath: output.decompositionData,
          err: output.decompositionError,
        });

        mergeMsBurden({
          plotPath: output.burdenPath,
          err: output.burdenError,
        });

        mergeMsAssociation({
          plotPath: output.associationPath,
          err: output.associationError,
        });

        mergeMsLandscape({
          plotPath: output.landscapePath,
          err: output.landscapeError,
        });

        mergeMsPrevalence({
          plotPath: output.prevalencePath,
          err: output.prevalenceError,
        });

        mergeMsIndividual({
          plotPath: output.individualPath,
          err: output.individualError,
        });
        mergeState({ submitted: true });
        if (displayTab == 'instructions')
          mergeState({ displayTab: 'tmb', openSidebar: false });
      } else {
        mergeError('');
      }
    } catch (err) {
      mergeError(err.message);
    }
    mergeState({ loading: false });
  }

  // when using public data and only need to upload a variable data file
  // also used to create a work directory and id
  async function uploadVariable() {
    return new Promise(async (resolve, reject) => {
      try {
        const data = new FormData();
        if (variableFileObj.size) data.append('variableFile', variableFileObj);

        let response = await fetch(`web/upload`, {
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
          await mergeState({ projectID });
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

  async function exposureDownload() {
    try {
      const { output, projectID, stdout } = await submitR('exposureDownload', {
        study: study,
        strategy: strategy,
        rsSet: rsSet,
        cancerType: cancer,
      });

      const file = await fetch(`web/results/${output.path}`);
      if (file.ok) {
        saveAs(await file.blob(), output.filename);
      } else {
        mergeError(`public data is not available`);
      }
    } catch (err) {
      console.log(err);
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

    const rsSetOptions = [
      ...new Set(
        exposureSignature
          .filter((row) => row.Study == study && row.Dataset == strategy)
          .map((row) => row.Signature_set_name)
      ),
    ];
    const rsSet = rsSetOptions[0];

    const cancerOptions = [
      ...new Set(
        exposureCancer
          .filter((data) => data.Study == study && data.Dataset == strategy)
          .map((data) => data.Cancer_Type)
      ),
    ];

    handleSet(rsSet);

    mergeState({
      study: study,
      strategy: strategy,
      strategyOptions: strategyOptions,
      cancer: cancerOptions[0],
      cancerOptions: cancerOptions,
      rsSetOptions: rsSetOptions,
      rsSet: rsSet,
    });
  }

  function handleStrategy(strategy) {
    const rsSetOptions = [
      ...new Set(
        exposureSignature
          .filter((row) => row.Study == study && row.Dataset == strategy)
          .map((row) => row.Signature_set_name)
      ),
    ];
    const rsSet = rsSetOptions[0];

    const cancerOptions = [
      ...new Set(
        exposureCancer
          .filter((data) => data.Study == study && data.Dataset == strategy)
          .map((data) => data.Cancer_Type)
      ),
    ];

    handleSet(rsSet);

    mergeState({
      strategy: strategy,
      cancer: cancerOptions[0],
      cancerOptions: cancerOptions,
      rsSetOptions: rsSetOptions,
      rsSet: rsSet,
    });
  }

  function handleSet(rsSet) {
    mergeState({ rsSet });
  }

  function handleReset() {
    window.location.hash = '#/exploration';

    const params = {
      source,
      study: 'PCAWG',
      strategy: 'WGS',
      rsSet: 'COSMIC_v3_Signatures_GRCh37_SBS96',
      cancer: 'Lung-AdenoCA',
      studyOptions,
      strategyOptions,
      rsSetOptions,
      cancerOptions,
      exposureCancer,
      exposureSignature,
    };
    resetExposure();
    mergeState(params);
  }

  const tabs =
    source == 'user'
      ? [
          {
            component: <Instructions loading={loading} />,
            id: 'instructions',
            name: 'Instructions',
          },
          {
            component: <TMB />,
            id: 'tmb',
            name: 'TMB',
          },
          {
            component: <TmbSig />,
            id: 'tmbSig',
            name: 'TMB Signatures',
          },
          {
            component: <MsBurden calculateBurden={calculateBurden} />,
            id: 'msBurden',
            name: 'MS Burden',
          },
          {
            component: <MsDecomposition />,
            id: 'msDecomposition',
            name: 'MS Decomposition',
          },
          {
            component: (
              <MsAssociation calculateAssociation={calculateAssociation} />
            ),
            id: 'msAssociation',
            name: 'MS Association',
          },
          {
            component: (
              <MsLandscape
                calculateLandscape={calculateLandscape}
                handleVariable={handleVariable}
              />
            ),
            id: 'msLandscape',
            name: 'MS Landscape',
          },
          {
            component: (
              <MsPrevalence calculatePrevalence={calculatePrevalence} />
            ),
            id: 'msPrevalence',
            name: 'MS Prevalence',
          },
          {
            component: (
              <MSIndividual calculateIndividual={calculateIndividual} />
            ),
            id: 'msIndividaul',
            name: 'MS Individual',
          },
          source == 'public' ? (
            {
              component: <Download exposureDownload={exposureDownload} />,
              id: 'download',
              name: 'Download',
            }
          ) : (
            <></>
          ),
        ]
      : [
          {
            component: <Instructions loading={loading} />,
            id: 'instructions',
            name: 'Instructions',
          },
          {
            component: <TMB2 />,
            id: 'tmb',
            name: 'TMB',
          },
          {
            component: <TmbSig2 />,
            id: 'tmbSig',
            name: 'TMB Signatures',
          },
          // {
          //   component: <MsBurden2 />,
          //   id: 'msBurden',
          //   name: 'MS Burden',
          // },
        ];

  return (
    <div className="position-relative">
      <div className="mx-3">
        <div className="mx-3 bg-white border border-top-0">
          {/* for desktops and tablets */}
          <div className="d-none d-md-block">
            <Nav defaultActiveKey="profilerSummary">
              {tabs.map(({ name, id }) => {
                if (name)
                  return (
                    <div key={id} className="d-inline-block">
                      <Button
                        variant="link"
                        className={`secondary-navlinks px-3 py-1 d-inline-block border-0 ${
                          id == displayTab ? 'active-secondary-navlinks' : ''
                        }`}
                        active={id == displayTab && submitted}
                        disabled={id != 'instructions' && !submitted}
                        style={{
                          textDecoration: 'none',
                          fontSize: '12pt',
                          color: 'black',
                          fontWeight: '500',
                        }}
                        onClick={() => mergeState({ displayTab: id })}
                      >
                        {name}
                      </Button>
                    </div>
                  );
              })}
            </Nav>
          </div>
          {/* for mobile devices */}
          <div className="row d-md-none">
            <Nav defaultActiveKey="summary">
              {tabs.map(({ name, id }) => {
                if (name)
                  return (
                    <div key={id} className="col-12 text-center">
                      <Button
                        variant="link"
                        className={
                          id == displayTab && Object.keys(exposureCancer).length
                            ? 'secondary-navlinks px-3 py-1 d-inline-block border-0 active-secondary-navlinks'
                            : 'secondary-navlinks px-3 py-1 d-inline-block border-0'
                        }
                        style={{
                          textDecoration: 'none',
                          fontSize: '12pt',
                          color: 'black',
                          fontWeight: '500',
                        }}
                        onClick={() => mergeState({ displayTab: id })}
                      >
                        {name}
                      </Button>
                      <div className="d-md-none w-100"></div>
                    </div>
                  );
              })}
            </Nav>
          </div>
        </div>
      </div>
      <SidebarContainer
        className="m-3"
        collapsed={!openSidebar}
        onCollapsed={(e) => mergeState({ openSidebar: !e })}
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
                      onChange={(e) => mergeState({ source: 'public' })}
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
                      onChange={(e) => mergeState({ source: 'user' })}
                    />
                    <Check.Label className="font-weight-normal">
                      User
                    </Check.Label>
                  </Check>
                </Group>
              </Col>
            </Row>
            <Row>
              <Col lg="12" className="w-100">
                {source == 'public' ? (
                  <PublicForm
                    calculate={handleCalculate}
                    handleReset={handleReset}
                    handleStudy={handleStudy}
                    handleStrategy={handleStrategy}
                    handleSet={handleSet}
                  />
                ) : (
                  <UserForm
                    calculate={handleCalculate}
                    handleReset={handleReset}
                    handleStudy={handleStudy}
                    handleStrategy={handleStrategy}
                  />
                )}
              </Col>
            </Row>
          </div>
        </SidebarPanel>
        <MainPanel>
          {tabs.filter((tab) => tab.id == displayTab)[0].component}
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
