import React, { useState, useEffect, useRef } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { actions as associationActions } from '../../../services/store/association';
import { actions as modalActions } from '../../../services/store/modal';
import Select from '../../controls/select/select';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

import Plot from '../../controls/plot/plot';
import Table from '../../controls/table/table';

const actions = { ...associationActions, ...modalActions };
const { Group, Label, Check, Control } = Form;

export default function Univariate() {
  const dispatch = useDispatch();
  const mergeState = async (state) =>
    await dispatch(actions.mergeAssociation({ univariate: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const {
    loadingData,
    assocVarData,
    expVarList,
    study,
    strategy,
    rsSet,
    cancer,
    assocTable,
  } = useSelector((state) => state.association.associationState);

  const {
    loadingParams,
    loadingCalculate,
    loadingRecalculate,
    error,
    assocVariant,
    expVariant,
    projectID,
    plotPath,
    dataPath,
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
    resultsTable,
  } = useSelector((state) => state.association.univariate);

  const [loadingMsg, setLoadingMsg] = useState(null);

  // populate controls
  useEffect(() => {
    if (assocVarData.length) {
      const dataSourceOptions = [
        ...new Set(assocVarData.map((row) => row.data_source)),
      ];
      const dataSource = dataSourceOptions[0];

      mergeState({
        dataSource,
        dataSourceOptions,
        expVariant: expVarList[0],
      });
    }
  }, [assocVarData]);
  // filter dropdown options on change
  useEffect(() => {
    if (dataSource) handleDataSource();
  }, [dataSource]);
  useEffect(() => {
    if (dataType) handleDataType();
  }, [dataType]);

  // recalculate on new signature name selection
  // ignore if prevSignature was '' (first calculation)
  const prevSignature = usePrevious(signature);
  useEffect(() => {
    if (signature && prevSignature != '') handleCalculate();
  }, [signature]);

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

  // hook for saving old values from previous render
  function usePrevious(value) {
    const ref = useRef();
    useEffect(() => {
      ref.current = value;
    });
    return ref.current;
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

  async function handleLoadParameters() {
    mergeState({ loadingParams: true, error: false });
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
        uncaught_error,
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

      if (uncaught_error) {
        console.error('R error: ' + uncaught_error);
        mergeState({
          error:
            'An error has occured. Please review your parameters and try again.',
        });
      } else if (error) {
        throw error.error;
      } else {
        mergeState({
          projectID: id,
          plotPath,
          dataPath,
          resultsTable: { data: dataTable },
          signatureOptions,
          signature: signature ? signature : signatureOptions[0],
        });
      }
    } catch (error) {
      mergeState({ error: error });
    }
    mergeState({
      loadingCalculate: false,
      loadingRecalculate: false,
    });
  }

  return (
    <div className="bg-white border rounded">
      <LoadingOverlay
        active={loadingData}
        content={loadingMsg}
        showIndicator={loadingMsg}
      />
      <div>
        <div className="mx-auto py-3 px-4">
          <Table
            title="Association Variable Data"
            data={assocVarData}
            columns={[
              ...new Set(...assocVarData.map((row) => Object.keys(row))),
            ].reduce(reducer, [])}
            pagination={assocTable.pagination}
            hidden={assocTable.hidden}
            mergeState={async (e) =>
              await dispatch(
                actions.mergeAssociation({ association: { assocTable: e } })
              )
            }
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
                    disabled={loadingData || loadingParams || loadingCalculate}
                    id="dataSource"
                    label="Data Source"
                    value={dataSource}
                    options={dataSourceOptions}
                    onChange={(e) => mergeState({ dataSource: e })}
                  />
                </Col>
                <Col md="4">
                  <Select
                    disabled={loadingData || loadingParams || loadingCalculate}
                    id="dataType"
                    label="Data Type"
                    value={dataType}
                    options={dataTypeOptions}
                    onChange={(e) => mergeState({ dataType: e })}
                  />
                </Col>
                <Col md="4">
                  <Select
                    disabled={loadingData || loadingParams || loadingCalculate}
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
                    disabled={loadingData || loadingParams || loadingCalculate}
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
                          loadingData || loadingParams || loadingCalculate
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
                              loadingData || loadingParams || loadingCalculate
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
                              loadingData || loadingParams || loadingCalculate
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
                              loadingData || loadingParams || loadingCalculate
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
                              loadingData || loadingParams || loadingCalculate
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
                    disabled={loadingData || loadingParams || loadingCalculate}
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
              Select an Association Variant and Signature Exposure Variant
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
                    ...resultsTable.data.map((row) => Object.keys(row))
                  ),
                ].reduce(reducer, [])}
                pagination={resultsTable.pagination}
                hidden={resultsTable.hidden}
                mergeState={async (e) =>
                  await mergeState({ resultsTable: { ...e } })
                }
              />
              <Row>
                <Col md="auto">
                  <p>
                    Select a Signature Name to recalculate for a different
                    Signature
                  </p>
                </Col>
                <Col md="auto">
                  <Select
                    disabled={loadingData || loadingParams || loadingCalculate}
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
                  height="800px"
                />
              )}
            </div>
          )}
        </div>
        <hr />
      </div>
    </div>
  );
}
