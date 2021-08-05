import React, { useState, useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
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

export default function Association() {
  const dispatch = useDispatch();
  const mergeState = async (state) =>
    await dispatch(actions.mergeAssociation({ ...state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));
  const resetAssociation = (_) => dispatch(actions.resetAssociation());

  const {
    openSidebar,
    loadingData,
    loadingParams,
    loadingCalculate,
    loadingRecalculate,
    submitted,
    error,
    exposureSignature,
    assocVarData,
    expVarList,
    assocVariant,
    expVariant,
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
    signature,
    signatureOptions,
    testType,
    xlab,
    ylab,
    variant1,
    variant2,
    assocTable,
    resultsTable,
  } = useSelector((state) => state.association);

  const [loadingMsg, setLoadingMsg] = useState(null);

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

  // recalculate on new signature name selection
  useEffect(() => {
    if (signature) handleCalculate();
  }, [signature]);

  // popualte side panel
  async function populateControls() {
    mergeState({ loadingData: true });

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
    } catch (error) {
      mergeError(error);
    }

    mergeState({ loadingData: false });
  }

  // reducer for creating table columns from objects
  const reducer = (acc, column) => [
    ...acc,
    {
      Header: column
        .replace('_', ' ')
        .replace(/(^\w{1})|(\s+\w{1})/g, (letter) => letter.toUpperCase()),
      accessor: column,
      id: column,
      Cell: (e) => e.value || '',
    },
  ];

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

    mergeState({ assocVariant: assocVarOptions[0], assocVarOptions });
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

    mergeState({ loadingData: true });
    try {
      const [assocVarData, exposureVariantData] = await Promise.all([
        getJSON(`Association/PCAWG_vardata.json`),
        loadData(),
      ]);

      const dataSourceOptions = [
        ...new Set(assocVarData.map((row) => row.data_source)),
      ];
      const dataSource = dataSourceOptions[0];

      const { expVarList, error } = exposureVariantData;
      if (error) throw error.message;
      
      mergeState({
        assocVarData,
        expVarList,
        dataSource,
        dataSourceOptions,
        expVariant: expVarList[0],
      });
    } catch (error) {
      mergeError(error);
    }
    mergeState({ loadingData: false });
  }

  async function handleLoadParameters() {
    mergeState({ loadingParams: true });
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
              assocVariant,
              expVariant,
            },
          }),
        })
      ).json();

      const { collapseVar1, collapseVar2 } = collapseData;

      mergeState({
        variant1: {
          name: assocVariant,
          collapseOptions: collapseVar1 || [],
        },
        variant2: {
          name: expVariant,
          collapseOptions: collapseVar2 || [],
        },
      });
    } catch (error) {
      mergeError(error);
    }
    mergeState({ loadingParams: false });
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
    mergeState({
      loadingCalculate: signature ? false : true,
      loadingRecalculate: signature ? true : false,
      error: false,
      plotPath: '',
      dataPath: '',
    });
    try {
      const {
        stdout,
        error,
        plotPath,
        dataPath,
        dataTable,
        signatureOptions,
        projectID: id,
      } = await (
        await fetch(`api/associationCalc`, {
          method: 'POST',
          headers: {
            Accept: 'application/json',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            fn: 'univariate',
            projectID,
            args: {
              study,
              strategy,
              rsSet,
              cancer,
              dataSource,
              dataType,
              testType,
              signature,
              xlab: xlab || variant1.name,
              ylab: ylab || variant2.name,
              variant1: (() => {
                const { collapseOptions, ...params } = variant1;
                return params;
              })(),
              variant2: (() => {
                const { collapseOptions, ...params } = variant2;
                return params;
              })(),
            },
          }),
        })
      ).json();

      if (error) {
        throw error.error;
      }
      mergeState({
        projectID: id,
        plotPath,
        dataPath,
        resultsTable: { data: dataTable },
        signatureOptions,
        signature: signature ? signature : signatureOptions[0],
      });
    } catch (error) {
      mergeState({ error: error });
    }
    mergeState({
      loadingCalculate: false,
      loadingRecalculate: false,
    });
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
                        disabled={
                          loadingData ||
                          loadingParams ||
                          loadingCalculate ||
                          projectID
                        }
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
                        disabled={
                          loadingData ||
                          loadingParams ||
                          loadingCalculate ||
                          projectID
                        }
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
                          disabled={
                            loadingData ||
                            loadingParams ||
                            loadingCalculate ||
                            submitted
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
                            loadingData ||
                            loadingParams ||
                            loadingCalculate ||
                            submitted
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
                            loadingData ||
                            loadingParams ||
                            loadingCalculate ||
                            submitted
                          }
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
                          disabled={
                            loadingData ||
                            loadingParams ||
                            loadingCalculate ||
                            submitted
                          }
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
                        disabled={
                          loadingData ||
                          loadingParams ||
                          loadingCalculate ||
                          loadingRecalculate
                        }
                        className="w-100 mb-3"
                        variant="secondary"
                        onClick={() => handleReset()}
                      >
                        Reset
                      </Button>
                    </Col>
                    <Col md="6">
                      <Button
                        disabled={
                          loadingData ||
                          loadingParams ||
                          loadingCalculate ||
                          submitted
                        }
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
                        disabled={loadingData || loadingParams || loadingCalculate || loadingRecalculate || submitted}
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
                        disabled={loadingData || loadingParams || loadingCalculate || loadingRecalculate || submitted}
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
                          disabled={loadingData || loadingParams || loadingCalculate || loadingRecalculate || submitted}
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
                            disabled={loadingData || loadingParams || loadingCalculate || loadingRecalculate || submitted}
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
                            disabled={loadingData || loadingParams || loadingCalculate || loadingRecalculate || submitted}
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
                          disabled={loadingData || loadingParams || loadingCalculate || loadingRecalculate || submitted}
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
                        disabled={loadingData || loadingParams || loadingCalculate || loadingRecalculate || submitted}
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
                      disabled={loadingData || loadingParams || loadingCalculate || loadingRecalculate }
                      className="w-100 mb-3"
                      variant="secondary"
                      onClick={() => handleReset()}
                    >
                      Reset
                    </Button>
                  </Col>
                  <Col md="6">
                    <Button
                      disabled={loadingData || loadingParams || loadingCalculate || loadingRecalculate }
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
                active={loadingData}
                content={loadingMsg}
                showIndicator={loadingMsg}
              />
              {assocVarData.length ? (
                <div>
                  <div className="mx-auto py-3 px-4">
                    <Table
                      title="Association Variable Data"
                      data={assocVarData}
                      columns={[
                        ...new Set(
                          ...assocVarData.map((row) => Object.keys(row))
                        ),
                      ].reduce(reducer, [])}
                      pagination={assocTable.pagination}
                      hidden={assocTable.hidden}
                      mergeState={(e) => mergeState({ assocTable: { ...e } })}
                    />
                  </div>
                  <hr />
                  <div className="mx-auto py-3 px-4">
                    <LoadingOverlay
                      active={loadingParams}
                      content={loadingMsg}
                      showIndicator={loadingMsg}
                    />
                    <h4>Select Variables</h4>
                    <Row className="justify-content-center mt-3">
                      <Col md="8">
                        <strong>Association Variable</strong>
                        <Row>
                          <Col md="4">
                            <Select
                              disabled={
                                loadingData ||
                                loadingParams ||
                                loadingCalculate ||
                                submitted
                              }
                              id="dataSource"
                              label="Data Source"
                              value={dataSource}
                              options={dataSourceOptions}
                              onChange={(e) => mergeState({ dataSource: e })}
                            />
                          </Col>
                          <Col md="4">
                            <Select
                              disabled={
                                loadingData ||
                                loadingParams ||
                                loadingCalculate ||
                                submitted
                              }
                              id="dataType"
                              label="Data Type"
                              value={dataType}
                              options={dataTypeOptions}
                              onChange={(e) => mergeState({ dataType: e })}
                            />
                          </Col>
                          <Col md="4">
                            <Select
                              disabled={
                                loadingData ||
                                loadingParams ||
                                loadingCalculate ||
                                submitted
                              }
                              id="assocVariant"
                              label="Variant Name"
                              value={assocVariant}
                              options={assocVarOptions}
                              onChange={(e) => mergeState({ assocVariant: e })}
                            />
                          </Col>
                        </Row>
                      </Col>
                      <Col md="4">
                        <strong>Signature Exposure Variable</strong>
                        <Row>
                          <Col md="12">
                            <Select
                              disabled={
                                loadingData ||
                                loadingParams ||
                                loadingCalculate ||
                                submitted
                              }
                              id="expVariant"
                              label="Variant Name"
                              value={expVariant}
                              options={expVarList}
                              onChange={(e) => mergeState({ expVariant: e })}
                            />
                          </Col>
                        </Row>
                      </Col>
                    </Row>
                    <Row className="justify-content-end">
                      <Col md="auto">
                        <Button
                          disabled={
                            loadingData ||
                            loadingParams ||
                            loadingCalculate ||
                            submitted ||
                            !dataSource
                          }
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
                    <LoadingOverlay
                      active={loadingCalculate}
                      content={loadingMsg}
                      showIndicator={loadingMsg}
                    />
                    <h4>Parameters</h4>
                    {variant1.name && variant2.name ? (
                      <>
                        <Row className="justify-content-center mt-3">
                          <Col md="12">
                            <Row>
                              <Col md="2">
                                <Select
                                  disabled={
                                    loadingData ||
                                    loadingParams ||
                                    loadingCalculate ||
                                    submitted
                                  }
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
                          <Col md="12">
                            <strong>Association Variant</strong>
                            <Row>
                              <Col md="auto">
                                <Group controlId="filter1" className="d-flex">
                                  <Label className="mr-4">Filtering (>0)</Label>
                                  <Check inline id="filter1">
                                    <Check.Input
                                      disabled={
                                        loadingData ||
                                        loadingParams ||
                                        loadingCalculate ||
                                        submitted
                                      }
                                      type="checkbox"
                                      value={variant1.filter}
                                      checked={variant1.filter}
                                      onChange={() =>
                                        mergeState({
                                          variant1: {
                                            filter: !variant1.filter,
                                          },
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
                                      disabled={
                                        loadingData ||
                                        loadingParams ||
                                        loadingCalculate ||
                                        submitted
                                      }
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
                                    loadingData ||
                                    loadingParams ||
                                    loadingCalculate ||
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
                          <Col md="12">
                            <strong>Signature Exposure Variant</strong>
                            <Row>
                              <Col md="auto">
                                <Group controlId="filter2" className="d-flex">
                                  <Label className="mr-4">Filtering (>0)</Label>
                                  <Check inline id="filter2">
                                    <Check.Input
                                      disabled={
                                        loadingData ||
                                        loadingParams ||
                                        loadingCalculate ||
                                        submitted
                                      }
                                      type="checkbox"
                                      value={variant2.filter}
                                      checked={variant2.filter}
                                      onChange={() =>
                                        mergeState({
                                          variant2: {
                                            filter: !variant2.filter,
                                          },
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
                                      disabled={
                                        loadingData ||
                                        loadingParams ||
                                        loadingCalculate ||
                                        submitted
                                      }
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
                                    loadingData ||
                                    loadingParams ||
                                    loadingCalculate ||
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
                              disabled={
                                loadingData ||
                                loadingParams ||
                                loadingCalculate ||
                                submitted
                              }
                              className="w-100"
                              variant="primary"
                              onClick={() => handleCalculate()}
                            >
                              Calculate
                            </Button>
                          </Col>
                        </Row>
                      </>
                    ) : (
                      <p className="d-flex justify-content-center text-muted">
                        Select an Association Variant and Signature Exposure
                        Variant
                      </p>
                    )}
                  </div>
                  <div>
                    <hr />
                    {error && (
                      <p className="p-3 d-flex justify-content-center text-danger">
                        {error}
                      </p>
                    )}
                    {resultsTable.data.length > 0 && (
                      <div className="mx-auto p-3">
                        <h4>Results</h4>
                        <Table
                          title="Association Group"
                          data={resultsTable.data}
                          columns={[
                            ...new Set(
                              ...resultsTable.data.map((row) =>
                                Object.keys(row)
                              )
                            ),
                          ].reduce(reducer, [])}
                          pagination={resultsTable.pagination}
                          hidden={resultsTable.hidden}
                          mergeState={async (e) =>
                            await mergeState({ resultsTable: { ...e } })
                          }
                        />
                        <Row>
                          <Col>
                            <p>
                              Select an Association Variant and Signature
                              Exposure Variant
                            </p>
                          </Col>
                          <Col md="auto">
                            <Select
                              disabled={
                                loadingData ||
                                loadingParams ||
                                loadingCalculate ||
                                submitted
                              }
                              id="signature"
                              label="Signature Name"
                              value={signature}
                              options={signatureOptions}
                              onChange={(e) => mergeState({ signature: e })}
                            />
                          </Col>
                        </Row>
                        <LoadingOverlay
                          active={loadingRecalculate}
                          content={loadingMsg}
                          showIndicator={loadingMsg}
                        />
                        {plotPath && (
                          <Plot
                            className="p-3 border rounded"
                            // title="Association"
                            downloadName={plotPath.split('/').slice(-1)[0]}
                            plotPath={`api/results/${projectID}${plotPath}`}
                            txtPath={projectID + dataPath}
                            maxHeight="800px"
                          />
                        )}
                      </div>
                    )}
                    {/* <Debug msg={debugR} /> */}
                  </div>
                  <hr />
                </div>
              ) : (
                <div className=" py-3 px-4">
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
    </div>
  );
}
