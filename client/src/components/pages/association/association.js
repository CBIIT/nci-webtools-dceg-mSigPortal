import React, { useState, useEffect, useRef } from 'react';
import { Form, Row, Col, Tab, Button, Nav } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import {
  actions as associationActions,
  getInitialState,
} from '../../../services/store/association';
import { actions as modalActions } from '../../../services/store/modal';
import Select from '../../controls/select/select';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '../../controls/sidebar-container/sidebar-container';
import Plot from '../../controls/plot/plot';
import Table from '../../controls/table/table';
import { getJSON } from '../../../services/utils';
import './association.scss';

const actions = { ...associationActions, ...modalActions };
const { Group, Label, Check, Control } = Form;

export default function Association({ match }) {
  const dispatch = useDispatch();
  const mergeState = async (state) =>
    await dispatch(actions.mergeAssociation({ ...state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));
  const resetAssociation = (_) => dispatch(actions.resetAssociation());

  const {
    openSidebar,
    loading,
    submitted,
    err,
    exposureSignature,
    assocVarData,
    expVarList,
    projectID,
    plotPath,
    dataPath,
    source,
    study,
    studyOptions,
    strategy,
    strategyOptions,
    rsSet,
    rsSetOptions,
    cancer,
    cancerOptions,
    dataSource,
    dataSourceOptions,
    dataType,
    dataTypeOptions,

    assocVarOptions,
    regression,
    testType,
    xlab,
    ylab,
    variant1,
    variant2,
  } = useSelector((state) => state.association);

  const [loadingMsg, setLoadingMsg] = useState(null);
  const [display, setDisplay] = useState('association');

  // populate controls on inital render
  useEffect(() => {
    if (!studyOptions.length) populateControls();
  }, []);

  // filter dropdown options on change
  useEffect(() => {
    if (dataSource) handleDataSource();
  }, [dataSource]);
  useEffect(() => {
    if (dataType) handleDataType();
  }, [dataType]);

  // popualte side panel
  async function populateControls() {
    mergeState({ loading: true });

    try {
      const exposureSignature = await getJSON(
        'Others/json/Exploring-Exposure.json'
      );

      const studyOptions = [
        ...new Set(exposureSignature.map((data) => data.Study)),
      ];
      const study = 'PCAWG'; // default

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
      const rsSet = 'COSMIC v3 Signatures (SBS)'; // default

      const cancerOptions = [
        ...new Set(
          exposureSignature
            .filter((data) => data.Study == study && data.Dataset == strategy)
            .map((data) => data.Cancer_Type)
        ),
      ];
      const cancer = 'Lung-AdenoCA'; // default

      mergeState({
        exposureSignature,
        study,
        studyOptions,
        strategy,
        strategyOptions,
        cancer,
        cancerOptions,
        rsSet,
        rsSetOptions,
      });
    } catch (err) {
      mergeError(err.message);
    }

    mergeState({ loading: false });
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
        exposureSignature
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
        exposureSignature
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
    const rsSetOptions = [
      ...new Set(
        exposureSignature
          .filter((row) => row.Signature_set_name == set)
          .map((row) => row.Signature_name)
      ),
    ];

    mergeState({
      rsSet: set,
      rsSetOptions: rsSetOptions,
    });
  }

  function handleDataSource() {
    const dataTypeOptions = [
      ...new Set(
        assocVarData
          .filter((row) => row.data_source == dataSource)
          .map((row) => row.data_type)
      ),
    ];

    mergeState({ dataType: dataTypeOptions[0], dataTypeOptions });
  }

  function handleDataType() {
    const assocVarOptions = [
      ...new Set(
        assocVarData
          .filter((row) => row.data_type == dataType)
          .map((row) => row.variable_name)
      ),
    ];

    mergeState({ variant1: { name: assocVarOptions[0] }, assocVarOptions });
  }

  async function handleLoadData() {
    const loadData = async () =>
      (
        await fetch(`api/associationData`, {
          method: 'POST',
          headers: {
            Accept: 'application/json',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            fn: 'loadData',
            args: { study, strategy, rsSet, cancer },
          }),
        })
      ).json();

    mergeState({ loading: true });
    try {
      const [assocVarData, exposureVariantData] = await Promise.all([
        getJSON(`Association/PCAWG_vardata.json`),
        loadData(),
      ]);

      const dataSourceOptions = [
        ...new Set(assocVarData.map((row) => row.data_source)),
      ];
      const dataSource = dataSourceOptions[0];

      const { expVarList } = exposureVariantData.output;

      mergeState({
        assocVarData,
        expVarList,
        dataSource,
        dataSourceOptions,
        variant2: {
          name: expVarList[0],
        },
      });
    } catch (err) {
      mergeError(err.message);
    }
    mergeState({ loading: false });
  }

  async function handleLoadParameters() {
    mergeState({ loading: true });
    try {
      const collapseData = await (
        await fetch(`api/associationData`, {
          method: 'POST',
          headers: {
            Accept: 'application/json',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            fn: 'loadCollapse',
            args: {
              study,
              strategy,
              rsSet,
              cancer,
              dataSource,
              dataType,
              assocVar: variant1.name,
              expVar: variant2.name,
            },
          }),
        })
      ).json();

      const { collapseVar1, collapseVar2 } = collapseData.output;

      mergeState({
        variant1: {
          collapseOptions: collapseVar1 || [],
        },
        variant2: {
          collapseOptions: collapseVar2 || [],
        },
      });
    } catch (err) {
      mergeError(err.message);
    }
    mergeState({ loading: false });
  }

  function handleReset() {
    const params = {
      source,
      exposureSignature,
      studyOptions,
      strategyOptions,
      cancerOptions,
      rsSetOptions,
      study,
      strategy,
      cancer,
      rsSet,
    };
    resetAssociation();
    mergeState(params);
  }

  async function handleCalculate() {
    mergeState({ loading: true, err: false });
    try {
      const { debugR, output, projectID: id } = await (
        await fetch(`api/associationCalc`, {
          method: 'POST',
          headers: {
            Accept: 'application/json',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            fn: 'calculate',
            projectID,
            args: {
              study,
              strategy,
              rsSet,
              cancer,
              dataSource,
              dataType,
              assocVar: variant1.name,
              expVar: variant2.name,
              collapse1: variant1.collapse || null,
              collapse2: variant2.collapse || null,
              filter1: variant1.filter,
              filter2: variant2.filter,
              log2_1: variant1.log2,
              log2_2: variant2.log2,
              regression,
              testType,
              xlab: xlab || variant1.name,
              ylab: ylab || variant2.name,
            },
          }),
        })
      ).json();

      mergeState({
        projectID: id,
        plotPath: output.plotPath,
        dataPath: output.dataPath,
      });
    } catch (err) {
      mergeError(err.message);
      mergeState({ err: true });
    }
    mergeState({ loading: false });
  }

  return (
    <div className="px-0 m-3 position-relative">
      <div className="mx-3">
        <SidebarContainer
          collapsed={!openSidebar}
          onCollapsed={(e) => mergeState({ openSidebar: !e })}
        >
          <SidebarPanel>
            <div className="p-3 bg-white border rounded">
              <Row>
                <Col md="auto">
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
              {source == 'public' ? (
                <div>
                  <Row>
                    <Col>
                      <Group>
                        <Select
                          disabled={loading || submitted}
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
                          disabled={loading || submitted}
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
                          disabled={loading || submitted}
                          id="expSetPublic"
                          label="Reference Signature Set"
                          value={rsSet}
                          options={rsSetOptions}
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
                          disabled={loading || submitted}
                          id="prevalenceCancerType"
                          label="Cancer Type"
                          value={cancer}
                          options={cancerOptions}
                          onChange={(cancer) =>
                            mergeState({
                              cancer: cancer,
                            })
                          }
                        />
                      </Group>
                    </Col>
                  </Row>
                  <Row>
                    <Col md="6">
                      <Button
                        disabled={loading}
                        className="w-100 mb-3"
                        variant="secondary"
                        onClick={() => handleReset()}
                      >
                        Reset
                      </Button>
                    </Col>
                    <Col md="6">
                      <Button
                        disabled={loading || submitted}
                        className="w-100"
                        variant="primary"
                        onClick={() => handleLoadData()}
                      >
                        Load Data
                      </Button>
                    </Col>
                  </Row>
                </div>
              ) : (
                <div>
                  {/* <Row>
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
                          if (e.target.files.length) {
                            setExposure(e.target.files[0]);
                            mergeState({
                              exposureFile: e.target.files[0].name,
                            });
                          }
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
                          if (e.target.files.length) {
                            setMatrix(e.target.files[0]);
                            mergeState({
                              matrixFile: e.target.files[0].name,
                            });
                          }
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
                            mergeState({
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
                            value={rsSet}
                            options={rsSetOptions}
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
                            if (e.target.files.length) {
                              setSignature(e.target.files[0]);
                              mergeState({
                                signatureFile: e.target.files[0].name,
                              });
                            }
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
                        onChange={(genome) => mergeState({ genome: genome })}
                      />
                    </Group>
                  </Col>
                </Row>
                <Row>
                  <Col md="6">
                    <Button
                      disabled={loading}
                      className="w-100 mb-3"
                      variant="secondary"
                      onClick={() => handleReset()}
                    >
                      Reset
                    </Button>
                  </Col>
                  <Col md="6">
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
                </Row> */}
                </div>
              )}
            </div>
          </SidebarPanel>
          <MainPanel>
            <div className="bg-white border rounded">
              <LoadingOverlay
                active={loading}
                content={loadingMsg}
                showIndicator={loadingMsg}
              />
              <div>
                <div className="mx-auto py-3 px-4">
                  <span>Select Variables</span>
                  <Row className="justify-content-center mt-3">
                    <Col as="fieldset" md="8" className="border rounded">
                      <legend>Association Variable</legend>
                      <Row>
                        <Col md="4">
                          <Select
                            disabled={loading || submitted}
                            id="dataSource"
                            label="Data Source"
                            value={dataSource}
                            options={dataSourceOptions}
                            onChange={(e) => mergeState({ dataSource: e })}
                          />
                        </Col>
                        <Col md="4">
                          <Select
                            disabled={loading || submitted}
                            id="dataType"
                            label="Data Type"
                            value={dataType}
                            options={dataTypeOptions}
                            onChange={(e) => mergeState({ dataType: e })}
                          />
                        </Col>
                        <Col md="4">
                          <Select
                            disabled={loading || submitted}
                            id="variantName"
                            label="Variant Name"
                            value={variant1.name}
                            options={assocVarOptions}
                            onChange={(e) =>
                              mergeState({ variant1: { name: e } })
                            }
                          />
                        </Col>
                      </Row>
                    </Col>
                    <Col as="fieldset" md="4" className="border rounded">
                      <legend>Signature Exposure Variable</legend>
                      <Row>
                        <Col md="12">
                          <Select
                            disabled={loading || submitted}
                            id="sigExpVar"
                            label="Variant Name"
                            value={variant2.name}
                            options={expVarList}
                            onChange={(e) =>
                              mergeState({ variant2: { name: e } })
                            }
                          />
                        </Col>
                      </Row>
                    </Col>
                  </Row>
                  <Row className="justify-content-end">
                    <Col md="auto">
                      <Button
                        disabled={loading || submitted || !dataSource}
                        className="w-100"
                        variant="primary"
                        onClick={() => handleLoadParameters()}
                      >
                        Load Parameters
                      </Button>
                    </Col>
                  </Row>
                </div>
                <hr />
                <div className="mx-auto py-3 px-4">
                  <strong>Parameters</strong>
                  <Row className="justify-content-center mt-3">
                    <Col as="fieldset" md="12" className="border rounded">
                      <legend>Plot</legend>
                      <Row>
                        <Col md="auto">
                          <Group controlId="regression" className="d-flex">
                            <Label className="mr-4">Regression</Label>
                            <Check inline id="regression">
                              <Check.Input
                                disabled={loading || submitted}
                                type="checkbox"
                                value={regression}
                                checked={regression}
                                onChange={() =>
                                  mergeState({
                                    regression: !regression,
                                  })
                                }
                              />
                            </Check>
                          </Group>
                        </Col>
                        <Col md="2">
                          <Select
                            disabled={loading || submitted}
                            id="testType"
                            label="Test Type"
                            value={testType}
                            options={['nonparametric', 'parametric']}
                            onChange={(e) => mergeState({ testType: e })}
                          />
                        </Col>
                        <Col md="auto">
                          <Group controlId="xlab">
                            <Label>X-label</Label>
                            <Control
                              value={xlab}
                              placeholder={variant1.name}
                              onChange={(e) =>
                                mergeState({
                                  xlab: e.target.value,
                                })
                              }
                              isInvalid={false}
                            />
                            <Form.Control.Feedback type="invalid">
                              Enter a valid label
                            </Form.Control.Feedback>
                          </Group>
                        </Col>
                        <Col md="auto">
                          <Group controlId="ylab">
                            <Label>Y-label</Label>
                            <Control
                              value={ylab}
                              placeholder={variant2.name}
                              onChange={(e) =>
                                mergeState({
                                  ylab: e.target.value,
                                })
                              }
                              isInvalid={false}
                            />
                            <Form.Control.Feedback type="invalid">
                              Enter a valid label
                            </Form.Control.Feedback>
                          </Group>
                        </Col>
                      </Row>
                    </Col>
                    <Col />
                  </Row>
                  <Row className="justify-content-center">
                    <Col as="fieldset" md="12" className="border rounded">
                      <legend>Association Variant</legend>
                      <Row>
                        <Col md="auto">
                          <Group controlId="filter1" className="d-flex">
                            <Label className="mr-4">Filtering (>0)</Label>
                            <Check inline id="filter1">
                              <Check.Input
                                disabled={loading || submitted}
                                type="checkbox"
                                value={variant1.filter}
                                checked={variant1.filter}
                                onChange={() =>
                                  mergeState({
                                    variant1: { filter: !variant1.filter },
                                  })
                                }
                              />
                            </Check>
                          </Group>
                        </Col>
                        <Col md="auto">
                          <Group controlId="log2-1" className="d-flex">
                            <Label className="mr-4">
                              log<sub>2</sub>
                            </Label>
                            <Check inline id="log2-1">
                              <Check.Input
                                disabled={loading || submitted}
                                type="checkbox"
                                value={variant1.log2}
                                checked={variant1.log2}
                                onChange={() =>
                                  mergeState({
                                    variant1: { log2: !variant1.log2 },
                                  })
                                }
                              />
                            </Check>
                          </Group>
                        </Col>
                        <Col md="3">
                          <Select
                            disabled={
                              loading ||
                              submitted ||
                              !variant1.collapseOptions.length
                            }
                            id="collapse1"
                            label="Collapse"
                            value={
                              variant1.collapseOptions.length
                                ? variant1.collapse
                                : 'None'
                            }
                            options={variant1.collapseOptions}
                            onChange={(e) =>
                              mergeState({ variant1: { collapse: e } })
                            }
                          />
                        </Col>
                      </Row>
                    </Col>
                  </Row>
                  <Row className="justify-content-center">
                    <Col as="fieldset" md="12" className="border rounded">
                      <legend>Signature Exposure Variant</legend>
                      <Row>
                        <Col md="auto">
                          <Group controlId="filter2" className="d-flex">
                            <Label className="mr-4">Filtering (>0)</Label>
                            <Check inline id="filter2">
                              <Check.Input
                                disabled={loading || submitted}
                                type="checkbox"
                                value={variant2.filter}
                                checked={variant2.filter}
                                onChange={() =>
                                  mergeState({
                                    variant2: { filter: !variant2.filter },
                                  })
                                }
                              />
                            </Check>
                          </Group>
                        </Col>
                        <Col md="auto">
                          <Group controlId="log2-2" className="d-flex">
                            <Label className="mr-4">
                              log<sub>2</sub>
                            </Label>
                            <Check inline id="log2-2">
                              <Check.Input
                                disabled={loading || submitted}
                                type="checkbox"
                                value={variant2.log2}
                                checked={variant2.log2}
                                onChange={() =>
                                  mergeState({
                                    variant2: { log2: !variant2.log2 },
                                  })
                                }
                              />
                            </Check>
                          </Group>
                        </Col>
                        <Col md="3">
                          <Select
                            disabled={
                              loading ||
                              submitted ||
                              !variant2.collapseOptions.length
                            }
                            id="collapse2"
                            label="Collapse"
                            value={
                              variant2.collapseOptions.length
                                ? variant2.collapse
                                : 'None'
                            }
                            options={variant2.collapseOptions}
                            onChange={(e) =>
                              mergeState({ variant2: { collapse: e } })
                            }
                          />
                        </Col>
                      </Row>
                    </Col>
                  </Row>
                  <Row className="justify-content-end">
                    <Col md="auto">
                      <Button
                        disabled={loading || submitted}
                        className="w-100"
                        variant="primary"
                        onClick={() => handleCalculate()}
                      >
                        Calculate
                      </Button>
                    </Col>
                  </Row>
                </div>
                <hr />
                <div className="mx-auto p-3">
                  <strong>Results</strong>
                  <div id="exposureAssociationPlot">
                    {err && (
                      <div>
                        <hr />
                        <p className="p-3 text-danger">{err}</p>
                      </div>
                    )}
                    {plotPath && (
                      <>
                        <hr />
                        <Plot
                          className="p-3"
                          // title="Association"
                          downloadName={plotPath.split('/').slice(-1)[0]}
                          plotPath={`api/results/${projectID}${plotPath}`}
                          txtPath={projectID + dataPath}
                          maxHeight="800px"
                        />
                      </>
                    )}
                    {/* <Debug msg={debugR} /> */}
                  </div>
                </div>
                <hr />
              </div>
            </div>
          </MainPanel>
        </SidebarContainer>
      </div>
    </div>
  );
}
