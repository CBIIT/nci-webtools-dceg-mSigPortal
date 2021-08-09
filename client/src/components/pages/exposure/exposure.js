import React, { useState, useEffect } from 'react';
import { Form, Row, Col, Button, Nav } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import PublicForm from './publicForm';
import UserForm from './userForm';
import TMB from './tmb';
import TmbSig from './tmbSignatures';
import MsBurden from './msBurden';
import MsAssociation from './msAssociation';
import MsDecomposition from './msDecomposition';
import MsLandscape from './msLandscape';
import MsPrevalence from './msPrevalence';
import MSIndividual from './msIndividual';
import Download from './download';
import { actions as exposureActions } from '../../../services/store/exposure';
import { actions as modalActions } from '../../../services/store/modal';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
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
    dispatch(actions.mergeExposure({ exposureState: state }));
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
    signatureNames,
    study,
    studyOptions,
    strategy,
    strategyOptions,
    cancer,
    cancerOptions,
    rsSet,
    rsSetOptions,
    useCancerType,
    signatureNameOptions,
    publicSampleOptions,
    genome,
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
  } = exposureStore.exposureState;

  const { loading: loadingMsBurden, ...burdenArgs } = exposureStore.msBurden;
  const {
    loading: loadingMsAssociation,
    ...associationArgs
  } = exposureStore.msAssociation;
  const {
    loading: loadingMsLandscape,
    ...landscapeArgs
  } = exposureStore.msLandscape;
  const {
    loading: loadingMsPrevalence,
    ...prevalenceArgs
  } = exposureStore.msPrevalence;
  const {
    loading: loadingMsIndividual,
    ...individualArgs
  } = exposureStore.msIndividual;

  const [expand, setExpand] = useState(false);

  // load example if available
  useEffect(() => {
    if (exampleName) loadExample(exampleName);
  }, [exampleName]);

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
      const { state } = await (
        await fetch(`api/getExposureExample/${id}`)
      ).json();

      mergeState(state);
    } catch (error) {
      mergeError(error.message);
    }
    mergeState({
      exposure: {
        loading: false,
      },
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

        await handleCalculate('burden', projectID || (await uploadVariable()));

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

        await handleCalculate(
          'association',
          projectID || (await uploadVariable())
        );

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

        await handleCalculate(
          'landscape',
          projectID || (await uploadVariable())
        );

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

        await handleCalculate(
          'prevalence',
          projectID || (await uploadVariable())
        );

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

        await handleCalculate(
          'individual',
          projectID || (await uploadVariable())
        );

        mergeMsIndividual({ loading: false });
      }
    } catch (error) {
      mergeError(error.message);
    }
  }

  async function submitR(fn, args, id = projectID) {
    try {
      const response = await fetch(`api/explorationCalc`, {
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
    mergeState({ loading: true });

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
    if (
      !loadingMsBurden &&
      !gettingSignatureNames &&
      (fn == 'all' || fn == 'burden')
    ) {
      args.burden = JSON.stringify({
        signatureName: burdenArgs.signatureName || signatureNameOptions[0],
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
        signatureName1:
          associationArgs.signatureName1 || signatureNameOptions[0],
        signatureName2:
          associationArgs.signatureName2 ||
          signatureNameOptions[1] ||
          signatureNameOptions[0],
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
        sample: individualArgs.sample || publicSampleOptions[0],
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
      const { debugR, output, errors, projectID: pID } = await submitR(
        rFn,
        args,
        id
      );

      if (output) {
        mergeState({ projectID: pID });

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

        mergeMsBurden({
          plotPath: output.burdenPath,
          err: errors.burdenError,
          debugR: debugR,
        });

        mergeMsAssociation({
          plotPath: output.associationPath,
          err: errors.associationError,
          debugR: debugR,
        });

        mergeMsLandscape({
          plotPath: output.landscapePath,
          err: errors.landscapeError,
          debugR: debugR,
        });

        mergeMsPrevalence({
          plotPath: output.prevalencePath,
          err: errors.prevalenceError,
          debugR: debugR,
        });

        mergeMsIndividual({
          plotPath: output.individualPath,
          err: errors.individualError,
          debugR: debugR,
        });
        mergeState({ submitted: true });
      } else {
        mergeError(debugR);
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
      const { output, projectID, debugR } = await submitR('exposureDownload', {
        study: study,
        strategy: strategy,
        rsSet: rsSet,
        cancerType: cancer,
      });

      const file = await fetch(`api/results/${projectID}${output.path}`);
      if (file.ok) {
        const objectURL = URL.createObjectURL(await file.blob());
        const tempLink = document.createElement('a');

        tempLink.href = `${objectURL}`;
        tempLink.setAttribute('download', output.filename);
        document.body.appendChild(tempLink);
        tempLink.click();
        document.body.removeChild(tempLink);
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

  function handleSet(set) {
    const signatureNameOptions = [
      ...new Set(
        signatureNames
          .filter((row) => row.Signature_set_name == set)
          .map((row) => row.Signature_name)
      ),
    ];

    mergeState({
      rsSet: set,
      signatureNameOptions: signatureNameOptions,
    });
  }

  function handleReset() {
    window.location.hash = '#/exposure';

    const params = {
      source: source,
      study: 'PCAWG',
      strategy: 'WGS',
      rsSet: 'COSMIC v3 Signatures (SBS)',
      cancer: 'Lung-AdenoCA',
      studyOptions: studyOptions,
      strategyOptions: strategyOptions,
      rsSetOptions: rsSetOptions,
      cancerOptions: cancerOptions,
    };
    resetExposure();
    mergeState(params);
  }

  const tabs = [
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
      component: <MsAssociation calculateAssociation={calculateAssociation} />,
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
      component: <MsPrevalence calculatePrevalence={calculatePrevalence} />,
      id: 'msPrevalence',
      name: 'MS Prevalence',
    },
    {
      component: <MSIndividual calculateIndividual={calculateIndividual} />,
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
  ];

  const queries = [
    {
      study: 'PCAWG',
      examples: [
        {
          name: 'Lung',
          title:
            'PCAWG/WGS/COSMIC v3 Signatures (SBS)/ Lung-AdenoCA; MSA SBS5 vs SBS40',
          path: 'pcawg-lungadenoca',
        },
        {
          name: 'Skin',
          title:
            'PCAWG/WGS/COSMIC v3 Signatures (SBS)/ Skin-Melanoma; MSA SBS7a vs SBS7b',
          path: 'pcawg-skinmelanoma',
        },
        {
          name: 'Breast',
          title:
            'PCAWG/WGS/COSMIC v3 Signatures (SBS)/ Breast-AdenoCA; MSA SBS3 vs SBS5',
          path: 'pcawg-breastadenoca',
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
      <div className="mx-3">
        <div className="mx-3 bg-white border border-top-0">
          {/* for desktops and tablets */}
          <div className="d-none d-md-block">
            <Nav defaultActiveKey="profilerSummary">
              {tabs.map(({ name, id }) => (
                <div key={id} className="d-inline-block">
                  <Button
                    variant="link"
                    className={`secondary-navlinks px-3 py-1 d-inline-block border-0 ${
                      id == displayTab && submitted
                        ? 'active-secondary-navlinks'
                        : ''
                    }`}
                    active={id == displayTab && submitted}
                    disabled={!submitted}
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
              ))}
            </Nav>
          </div>
          {/* for mobile devices */}
          <div className="row d-md-none">
            <Nav defaultActiveKey="summary">
              {tabs.map(({ name, id }) => (
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
              ))}
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
            <strong>Example Queries</strong>
            {queries.map(({ study, examples }, index) => {
              return index == 0 ? (
                <Row key={`${study}-${index}`} className="mb-2">
                  <Col md="3">{study}:</Col>
                  {examples.map(({ name, title, path }) => (
                    <Col md="3" key={name + index}>
                      <span className="mb-2" disabled={loading || submitted}>
                        <a href={`#/exposure/${path}`} title={title}>
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
          <div className="bg-white border rounded">
            {submitted ? (
              tabs.filter((tab) => tab.id == displayTab)[0].component
            ) : (
              <div className="py-3 px-4">
                <LoadingOverlay active={loading} />
                <h4>Instructions</h4>
                <p>
                  Choose a Data Source and its associated options to submit a
                  query using the panel on the left
                </p>
                <hr />
                <h4>Data Source</h4>
                <p>
                  Public: Perform analysis using data available on the website
                </p>
                <p>User: Upload your own data</p>
              </div>
            )}
          </div>
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
