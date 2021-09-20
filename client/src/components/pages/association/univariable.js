import React, { useEffect, useRef, useMemo } from 'react';
import {
  Form,
  Row,
  Col,
  Button,
  OverlayTrigger,
  Popover,
} from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faInfoCircle } from '@fortawesome/free-solid-svg-icons';
import { actions as associationActions } from '../../../services/store/association';
import { actions as modalActions } from '../../../services/store/modal';
import Select from '../../controls/select/select';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

import Plot from '../../controls/plot/plot';
import Table from '../../controls/table/table';

const actions = { ...associationActions, ...modalActions };
const { Group, Label, Check, Control } = Form;

export default function Univariable() {
  const dispatch = useDispatch();
  const mergeState = async (state) =>
    await dispatch(actions.mergeAssociation({ univariable: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const {
    loadingData,
    assocVarData,
    assocFullDataPath,
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
    assocVariable,
    expVariable,
    projectID,
    plotPath,
    dataPath,
    assocTablePath,
    variableSource,
    variableSourceOptions,
    variableType,
    variableTypeOptions,
    assocVarOptions,
    signature,
    signatureOptions,
    testType,
    xlab,
    ylab,
    variable1,
    variable2,
    resultsTable,
  } = useSelector((state) => state.association.univariable);

  // populate controls
  useEffect(() => {
    if (assocVarData.length && !variableSource) {
      const variableSourceOptions = [
        ...new Set(assocVarData.map((row) => row.data_source)),
      ];
      const variableSource = variableSourceOptions[0];

      mergeState({
        variableSource,
        variableSourceOptions,
        expVariable: expVarList[0],
      });
    }
  }, [assocVarData]);
  // filter dropdown options on change
  useEffect(() => {
    if (variableSource) handleVariableSource();
  }, [variableSource]);
  useEffect(() => {
    if (variableType) handleVariableType();
  }, [variableType]);

  // recalculate on new signature name selection
  // ignore if prevSignature was '' (first calculation)
  const prevSignature = usePrevious(signature);
  useEffect(() => {
    if (signature && prevSignature && signature != prevSignature)
      handleCalculate();
  }, [signature]);

  // reducer for creating table columns from objects
  const reducer = (acc, column) => [
    ...acc,
    {
      Header: column,
      id: column,
      accessor: (a) => a[column],
      // sortType: useMemo(() => (rowA, rowB, columnId) => {
      //   const a = Number(rowA.original[columnId]);
      //   const b = Number(rowB.original[columnId]);
      //   if (a > b) return 1;
      //   if (b > a) return -1;
      //   return 0;
      // }),
      // sortMethod: (a, b) => Number(a) - Number(b),
      // Cell: (e) => e.value || '',
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

  function handleVariableSource() {
    const variableTypeOptions = [
      ...new Set(
        assocVarData
          .filter((row) => row.data_source == variableSource)
          .map((row) => row.data_type)
      ),
    ];

    mergeState({ variableType: variableTypeOptions[0], variableTypeOptions });
  }

  function handleVariableType() {
    const assocVarOptions = [
      ...new Set(
        assocVarData
          .filter((row) => row.data_type == variableType)
          .map((row) => row.variable_name)
      ),
    ];

    mergeState({ assocVariable: assocVarOptions[0], assocVarOptions });
  }

  async function handleLoadParameters() {
    mergeState({ loadingParams: true, error: false });
    try {
      const collapseData = await (
        await fetch(`api/associationWrapper`, {
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
              variableSource,
              variableType,
              assocVariable,
              expVariable,
            },
          }),
        })
      ).json();

      const { collapseVar1, collapseVar2 } = collapseData;

      mergeState({
        variable1: {
          name: assocVariable,
          collapseOptions: collapseVar1 || [],
        },
        variable2: {
          name: expVariable,
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
        assocTablePath,
        dataTable,
        signatureOptions,
        projectID: id,
      } = await (
        await fetch(`api/associationWrapper`, {
          method: 'POST',
          headers: {
            Accept: 'application/json',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            fn: 'univariable',
            projectID,
            args: {
              study,
              strategy,
              rsSet,
              cancer,
              variableSource,
              variableType,
              testType,
              signature,
              xlab: xlab || variable1.name,
              ylab: ylab || variable2.name,
              variable1: (() => {
                const { collapseOptions, ...params } = variable1;
                return params;
              })(),
              variable2: (() => {
                const { collapseOptions, ...params } = variable2;
                return params;
              })(),
            },
          }),
        })
      ).json();

      if (error) {
        mergeState({ error: error.error });
      } else if (uncaught_error) {
        console.error('R error: ' + uncaught_error);
        mergeState({
          error:
            'An error has occured. Please review your parameters and try again.',
        });
      } else {
        mergeState({
          projectID: id,
          plotPath,
          dataPath,
          assocTablePath,
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

  const popoverInfo = (text) => (
    <Popover>
      <Popover.Content>{text}</Popover.Content>
    </Popover>
  );

  return (
    <div className="p-4 bg-white border rounded">
      <LoadingOverlay active={loadingData} />
      <div className="mb-3">
        <h5 className="separator">Variables</h5>
        <Table
          customTitle={
            <>
              <Col />
              <Col />
              <Col>Available variables for selected study</Col>
            </>
          }
          data={assocVarData}
          columns={[
            ...new Set(...assocVarData.map((row) => Object.keys(row))),
          ].reduce(reducer, [])}
          pagination={assocTable.pagination}
          hidden={assocTable.hidden}
          downloadName="Download Variable Data"
          downloadLink={projectID + assocFullDataPath}
          mergeState={async (e) =>
            await dispatch(
              actions.mergeAssociation({ association: { assocTable: e } })
            )
          }
        />
        <LoadingOverlay active={loadingParams} />
        <p className="text-center">
          Select the following variables for analysis
        </p>
        <Row className="justify-content-center mt-3">
          <Col md="auto">
            <Select
              disabled={loadingData || loadingParams || loadingCalculate}
              id="variableSource"
              label="Variable Source"
              value={variableSource}
              options={variableSourceOptions}
              onChange={(e) => mergeState({ variableSource: e })}
            />
          </Col>
          <Col md="auto" lg="3" xl="2">
            <Select
              disabled={loadingData || loadingParams || loadingCalculate}
              id="variableType"
              label="Data Type"
              value={variableType}
              options={variableTypeOptions}
              onChange={(e) => mergeState({ variableType: e })}
            />
          </Col>
          <Col md="auto" lg="3" xl="2">
            <Select
              disabled={loadingData || loadingParams || loadingCalculate}
              id="assocVariable"
              label="Variant Name"
              value={assocVariable}
              options={assocVarOptions}
              onChange={(e) => mergeState({ assocVariable: e })}
            />
          </Col>
          <Col lg="auto">
            <fieldset className="border rounded p-2">
              <legend className="font-weight-bold">Variable Filtering</legend>
              <Row>
                <Col md="auto">
                  <div className="d-flex">
                    <OverlayTrigger
                      trigger="click"
                      placement="top"
                      overlay={popoverInfo(
                        'Filter sample with variable value above this threshold'
                      )}
                      rootClose
                    >
                      <Button
                        aria-label="threshold info"
                        variant="link"
                        className="p-0 font-weight-bold mr-1"
                      >
                        <FontAwesomeIcon
                          icon={faInfoCircle}
                          style={{ verticalAlign: 'super' }}
                        />
                      </Button>
                    </OverlayTrigger>
                    <Group
                      controlId="variable1-threshold"
                      className="d-flex mb-0"
                    >
                      <Label className="mr-2">Threshold</Label>
                      <Control
                        disabled={
                          loadingData || loadingParams || loadingCalculate
                        }
                        value={variable1.filter}
                        placeholder={'optional'}
                        style={{ width: '90px' }}
                        onChange={(e) =>
                          mergeState({
                            variable1: {
                              filter: e.target.value,
                            },
                          })
                        }
                        isInvalid={false}
                      />
                      <Form.Control.Feedback type="invalid">
                        Enter a valid threshold
                      </Form.Control.Feedback>
                    </Group>
                  </div>
                </Col>
                <Col md="auto">
                  <Group controlId="log2-1" className="d-flex mb-0">
                    <OverlayTrigger
                      trigger="click"
                      placement="top"
                      overlay={popoverInfo([
                        'Log',
                        <sub>2</sub>,
                        ' transform variable value',
                      ])}
                      rootClose
                    >
                      <Button
                        aria-label="log info"
                        variant="link"
                        className="p-0 font-weight-bold mr-1"
                      >
                        <FontAwesomeIcon
                          icon={faInfoCircle}
                          style={{ verticalAlign: 'baseline' }}
                        />
                      </Button>
                    </OverlayTrigger>{' '}
                    <Label className="mr-2 font-weight-normal">
                      log<sub>2</sub>
                    </Label>
                    <Check inline id="log2-1">
                      <Check.Input
                        disabled={
                          loadingData || loadingParams || loadingCalculate
                        }
                        type="checkbox"
                        value={variable1.log2}
                        checked={variable1.log2}
                        onChange={() =>
                          mergeState({
                            variable1: { log2: !variable1.log2 },
                          })
                        }
                      />
                    </Check>
                  </Group>
                </Col>
                <Col md="auto">
                  <OverlayTrigger
                    trigger="click"
                    placement="top"
                    overlay={popoverInfo(
                      'If variable value is factor, group variable value in to select collapse level and other'
                    )}
                    rootClose
                  >
                    <Button
                      aria-label="collapse info"
                      variant="link"
                      className="p-0 font-weight-bold"
                    >
                      <FontAwesomeIcon
                        icon={faInfoCircle}
                        style={{ verticalAlign: 'baseline' }}
                      />
                    </Button>
                  </OverlayTrigger>{' '}
                  <Select
                    className="d-inline-flex mb-0"
                    disabled={
                      loadingData ||
                      loadingParams ||
                      loadingCalculate ||
                      !variable1.collapseOptions.length
                    }
                    id="collapse1"
                    label="Collapse"
                    labelClass="mr-2 font-weight-normal"
                    value={
                      variable1.collapseOptions.length
                        ? variable1.collapse
                        : 'None'
                    }
                    options={variable1.collapseOptions}
                    onChange={(e) => mergeState({ variable1: { collapse: e } })}
                  />
                </Col>
              </Row>
            </fieldset>
          </Col>
        </Row>
        <Row className="justify-content-end">
          <Col md="auto">
            <Button
              disabled={
                loadingData ||
                loadingParams ||
                loadingCalculate ||
                !variableSource
              }
              className="w-100"
              variant="primary"
              onClick={() => handleLoadParameters()}
            >
              Load Data
            </Button>
          </Col>
        </Row>
      </div>
      <div className="mb-3">
        <h5 className="separator">Parameters</h5>

        <LoadingOverlay active={loadingCalculate} />
        {variable1.name && variable2.name ? (
          <>
            <p className="text-center">
              Select the following filtering and method for analysis
            </p>
            <Row className="justify-content-center">
              <Col md="auto" lg="auto">
                <Select
                  disabled={loadingData || loadingParams || loadingCalculate}
                  id="expVariable"
                  label="Signature Exposure Variable"
                  value={expVariable}
                  options={expVarList}
                  onChange={(e) => mergeState({ expVariable: e })}
                />
              </Col>
              <Col lg="auto">
                <fieldset className="border rounded p-2">
                  <legend className="font-weight-bold">
                    Signature Exposure Filtering
                  </legend>
                  <Row>
                    <Col md="auto">
                      <div className="d-flex">
                        <OverlayTrigger
                          trigger="click"
                          placement="top"
                          overlay={popoverInfo(
                            'Filter sample with signature exposure value above this threshold'
                          )}
                          rootClose
                        >
                          <Button
                            aria-label="threshold info"
                            variant="link"
                            className="p-0 font-weight-bold mr-1"
                          >
                            <FontAwesomeIcon
                              icon={faInfoCircle}
                              style={{ verticalAlign: 'super' }}
                            />
                          </Button>
                        </OverlayTrigger>
                        <Group
                          controlId="variable2-threshold"
                          className="d-flex mb-0"
                        >
                          <Label className="mr-2">Threshold</Label>
                          <Control
                            disabled={
                              loadingData || loadingParams || loadingCalculate
                            }
                            value={variable2.filter}
                            placeholder={'optional'}
                            style={{ width: '90px' }}
                            onChange={(e) =>
                              mergeState({
                                variable2: {
                                  filter: e.target.value,
                                },
                              })
                            }
                            isInvalid={false}
                          />
                          <Form.Control.Feedback type="invalid">
                            Enter a valid threshold
                          </Form.Control.Feedback>
                        </Group>
                      </div>
                    </Col>
                    <Col md="auto">
                      <Group controlId="log2-2" className="d-flex mb-0">
                        <OverlayTrigger
                          trigger="click"
                          placement="top"
                          overlay={popoverInfo([
                            'Log ',
                            <sub>2</sub>,
                            ' transform signature exposure value',
                          ])}
                          rootClose
                        >
                          <Button
                            aria-label="log info"
                            variant="link"
                            className="p-0 font-weight-bold mr-1"
                          >
                            <FontAwesomeIcon
                              icon={faInfoCircle}
                              style={{ verticalAlign: 'baseline' }}
                            />
                          </Button>
                        </OverlayTrigger>
                        <Label className="mr-2 font-weight-normal">
                          log<sub>2</sub>
                        </Label>
                        <Check inline id="log2-2">
                          <Check.Input
                            disabled={
                              loadingData || loadingParams || loadingCalculate
                            }
                            type="checkbox"
                            value={variable2.log2}
                            checked={variable2.log2}
                            onChange={() =>
                              mergeState({
                                variable2: { log2: !variable2.log2 },
                              })
                            }
                          />
                        </Check>
                      </Group>
                    </Col>
                    {/* <Col md="3">
                      <Select
                        disabled={
                          loadingData ||
                          loadingParams ||
                          loadingCalculate ||
                          !variable2.collapseOptions.length
                        }
                        id="collapse2"
                        label="Collapse"
                        value={
                          variable2.collapseOptions.length
                            ? variable2.collapse
                            : 'None'
                        }
                        options={variable2.collapseOptions}
                        onChange={(e) =>
                          mergeState({ variable2: { collapse: e } })
                        }
                      />
                    </Col> */}
                  </Row>
                </fieldset>
              </Col>
              <Col lg="auto">
                <fieldset className="border rounded p-2">
                  <legend className="font-weight-bold">Method</legend>
                  <Select
                    className="mb-0"
                    disabled={loadingData || loadingParams || loadingCalculate}
                    id="testType"
                    label=""
                    value={testType}
                    options={['nonparametric', 'parametric']}
                    onChange={(e) => mergeState({ testType: e })}
                  />
                </fieldset>
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
            Select an Association Variable and Signature Exposure Variable
          </p>
        )}
      </div>
      {(resultsTable.data.length > 0 || error) && (
        <div className="mb-3">
          <h5 className="separator">Results</h5>
          {error && (
            <p className="p-3 d-flex justify-content-center text-danger">
              {error}
            </p>
          )}
          {resultsTable.data.length > 0 && (
            <>
              <div className="mb-4">
                <p className="text-center">
                  Check the following table for the analyses result for all
                  signatures
                </p>
                <Table
                  title=""
                  data={resultsTable.data}
                  columns={[
                    ...new Set(
                      ...resultsTable.data.map((row) => Object.keys(row))
                    ),
                  ].reduce(reducer, [])}
                  pagination={resultsTable.pagination}
                  hidden={resultsTable.hidden}
                  downloadName="Download Association Result"
                  downloadLink={projectID + assocTablePath}
                  mergeState={async (e) =>
                    await mergeState({ resultsTable: { ...e } })
                  }
                />
              </div>
              <p className="text-center">
                Select a Signature Name to recalculate for a different Signature
              </p>
              <Row className="justify-content-center">
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
                <Col md="auto">
                  <Group controlId="xlab">
                    <Label>Variable Title</Label>
                    <Control
                      value={xlab}
                      placeholder={variable1.name}
                      onChange={(e) => mergeState({ xlab: e.target.value })}
                      isInvalid={false}
                    />
                    <Form.Control.Feedback type="invalid">
                      Enter a valid label
                    </Form.Control.Feedback>
                  </Group>
                </Col>
                <Col md="auto">
                  <Group controlId="ylab">
                    <Label>Signature Exposure Title</Label>
                    <Control
                      value={ylab}
                      placeholder={variable2.name}
                      onChange={(e) => mergeState({ ylab: e.target.value })}
                      isInvalid={false}
                    />
                    <Form.Control.Feedback type="invalid">
                      Enter a valid label
                    </Form.Control.Feedback>
                  </Group>
                </Col>
              </Row>
              <LoadingOverlay active={loadingRecalculate} />
              {plotPath && (
                <Plot
                  className="p-3 border rounded"
                  downloadName={plotPath.split('/').slice(-1)[0]}
                  plotPath={`api/results/${projectID}${plotPath}`}
                  txtPath={projectID + dataPath}
                  height="800px"
                />
              )}
            </>
          )}
        </div>
      )}
    </div>
  );
}
