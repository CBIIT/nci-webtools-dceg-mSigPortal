import React, { useEffect, useRef, useState } from 'react';
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
import { faInfoCircle, faPlus } from '@fortawesome/free-solid-svg-icons';
import { actions as associationActions } from '../../../services/store/association';
import { actions as modalActions } from '../../../services/store/modal';
import Select from '../../controls/select/select';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import AssocVarParams from './assocVarParams';
import Plot from '../../controls/plot/plot';
import Table from '../../controls/table/table';

const actions = { ...associationActions, ...modalActions };
const { Group, Label, Check, Control } = Form;

export default function Multivariable() {
  const dispatch = useDispatch();
  const mergeState = async (state) =>
    await dispatch(actions.mergeAssociation({ multivariable: state }));
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

  const multivariable = useSelector((state) => state.association.multivariable);
  const {
    loadingParams,
    loadingCalculate,
    loadingRecalculate,
    error,
    projectID,
    plotPath,
    dataPath,
    assocTablePath,
    signature,
    signatureOptions,
    testType,
    xlab,
    ylab,
    associationVars,
    exposureVar,
    resultsTable,
  } = multivariable;

  const [warnLimit, showWarnLimit] = useState(false);
  const [warnDupe, setDupe] = useState([]);

  // populate controls
  useEffect(() => {
    if (expVarList.length && !exposureVar.name) {
      mergeState({ exposureVar: { name: expVarList[0] } });
    }
  }, [expVarList]);
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

  // returns a mapping of input params and an array of indexes for duplicate inputs
  function findDupes(assocVars = associationVars) {
    const params = assocVars.map(({ source, type, tmpName }) => ({
      source,
      type,
      name: tmpName,
    }));

    return params.reduce((a, e, i) => {
      const input = `${e.source}|${e.type}|${e.name}`;
      a[input] = a[input] ? [...a[input], i] : [];
      return a;
    }, {});
  }

  async function handleLoadData() {
    const dupeIndexes = Object.values(findDupes()).flat();
    setDupe(dupeIndexes);

    if (!dupeIndexes.length) {
      mergeState({ loadingParams: true, error: false });
      try {
        const { stdout, output: collapseData } = await (
          await fetch(`api/associationWrapper`, {
            method: 'POST',
            headers: {
              Accept: 'application/json',
              'Content-Type': 'application/json',
            },
            body: JSON.stringify({
              fn: 'loadCollapseMulti',
              args: {
                study,
                strategy,
                rsSet,
                cancer,
                expName: exposureVar.name,
                associationVars: associationVars.map(
                  ({
                    sourceOptions,
                    typeOptions,
                    nameOptions,
                    tmpName,
                    collapseOptions,
                    ...params
                  }) => ({ ...params, name: tmpName })
                ),
              },
            }),
          })
        ).json();

        mergeState({
          associationVars: associationVars.map((assocVar, i) => ({
            ...assocVar,
            name: assocVar.tmpName,
            collapse: Array.isArray(collapseData[i + 1])
              ? collapseData[i + 1][0]
              : '',
            collapseOptions: Array.isArray(collapseData[i + 1])
              ? collapseData[i + 1]
              : [],
          })),
          // exposureVar: { name: expVarList[0] },
        });
      } catch (error) {
        mergeError(error);
      }
      mergeState({ loadingParams: false });
    }
  }

  function handleReset() {
    mergeState({
      associationVars: [
        {
          source: '',
          type: '',
          tmpName: '',
          sourceOptions: [],
          typeOptions: [],
          nameOptions: [],
          filter: '',
          log2: false,
          collapse: '',
          collapseOptions: [],
        },
      ],
      projectID: '',
      plotPath: '',
      dataPath: '',
      assocTablePath: '',
      resultsTable: { data: [] },
      signatureOptions: [],
      signature: '',
      xlab: '',
      ylab: '',
      error: '',
    });
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
        projectID: id,
        stdout,
        output,
      } = await (
        await fetch(`api/associationWrapper`, {
          method: 'POST',
          headers: {
            Accept: 'application/json',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            fn: 'multivariable',
            projectID,
            args: {
              study,
              strategy,
              rsSet,
              cancer,
              testType,
              signature,
              xlab: xlab || associationVars.name,
              ylab: ylab || exposureVar.name,
              associationVars: associationVars.map(
                ({
                  sourceOptions,
                  typeOptions,
                  nameOptions,
                  tmpName,
                  collapseOptions,
                  ...params
                }) => params
              ),
              exposureVar: (() => {
                const { nameOptions, ...params } = exposureVar;
                return params;
              })(),
            },
          }),
        })
      ).json();

      const {
        plotPath,
        dataPath,
        assocTablePath,
        dataTable,
        signatureOptions,
        error,
        uncaughtError,
      } = output;

      if (error || uncaughtError) {
        mergeState({ error: error || uncaughtError });
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

  function addParam() {
    if (associationVars.length < 10) {
      let newParams = associationVars.slice();
      newParams.push({
        source: '',
        type: '',
        tmpName: '',
        sourceOptions: [],
        typeOptions: [],
        nameOptions: [],
        filter: '',
        log2: false,
        collapse: '',
        collapseOptions: [],
      });
      mergeState({ associationVars: newParams });
    } else {
      showWarnLimit(true);
      setTimeout(() => showWarnLimit(false), 5000);
    }
  }

  function removeParam(index) {
    let newParams = associationVars.slice();
    newParams.splice(index, 1);
    const dupeIndexes = Object.values(findDupes(newParams)).flat();
    setDupe(dupeIndexes);
    mergeState({ associationVars: newParams });
  }

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
          mergeState={(e) =>
            dispatch(
              actions.mergeAssociation({ association: { assocTable: e } })
            )
          }
        />
        <div>
          <LoadingOverlay active={loadingParams} />
          <p className="text-center">
            Select the following variables for analysis
          </p>

          {associationVars.map((paramState, index) => (
            <AssocVarParams
              index={index}
              hostState={multivariable}
              paramState={paramState}
              mergeState={(e) => {
                let newParams = associationVars.slice();
                newParams[index] = { ...newParams[index], ...e };
                mergeState({ associationVars: newParams });
              }}
              remove={index != 0 ? () => removeParam(index) : false}
              duplicates={warnDupe}
            />
          ))}

          <Row
            className="mx-auto mt-3 justify-content-between"
            style={{ maxWidth: '1720px' }}
          >
            <Col md="auto" className="d-flex">
              <OverlayTrigger
                show={warnLimit}
                placement="bottom"
                overlay={
                  <Popover>
                    <Popover.Content className="text-danger">
                      You may only add up to 10 variables
                    </Popover.Content>
                  </Popover>
                }
              >
                <Button
                  variant="link"
                  onClick={() => addParam()}
                  title="Add Plot"
                  style={{ textDecoration: 'none' }}
                >
                  <span
                    className="text-nowrap"
                    title="Add Association Variable"
                  >
                    <FontAwesomeIcon icon={faPlus} /> Add Association Variable
                  </span>
                </Button>
              </OverlayTrigger>
            </Col>
            <Col md="auto">
              {warnDupe.length > 0 && (
                <span className="text-danger">
                  Please change or remove duplicate variables
                </span>
              )}
            </Col>
            <Col md="auto">
              <Button
                disabled={loadingData}
                className="mr-4 reset"
                variant="secondary"
                onClick={() => handleReset()}
              >
                Reset
              </Button>
              <Button
                disabled={
                  loadingData ||
                  loadingParams ||
                  loadingCalculate ||
                  resultsTable.data.length
                }
                variant="primary"
                style={{ width: 'fit-content' }}
                onClick={() => handleLoadData()}
              >
                Load Data
              </Button>
            </Col>
          </Row>
        </div>
      </div>
      <div className="mb-3">
        <h5 className="separator">Parameters</h5>
        <LoadingOverlay active={loadingCalculate} />
        {associationVars.filter(({ name }) => name).length ==
          associationVars.length && exposureVar.name ? (
          <>
            <p className="text-center">
              Select the following filtering and method for analysis
            </p>
            <Row className="justify-content-center">
              <Col md="auto" lg="auto">
                <Select
                  disabled={
                    loadingData ||
                    loadingParams ||
                    loadingCalculate ||
                    resultsTable.data.length
                  }
                  id="expVariable"
                  label="Signature Exposure Variable"
                  value={exposureVar.name}
                  options={expVarList}
                  onChange={(e) => mergeState({ exposureVar: { name: e } })}
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
                          controlId="exposureVar-threshold"
                          className="d-flex mb-0"
                        >
                          <Label className="mr-2 font-weight-normal">
                            Threshold
                          </Label>
                          <Control
                            disabled={
                              loadingData ||
                              loadingParams ||
                              loadingCalculate ||
                              resultsTable.data.length
                            }
                            value={exposureVar.filter}
                            placeholder={'Optional'}
                            style={{ width: '90px' }}
                            onChange={(e) =>
                              mergeState({
                                exposureVar: {
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
                              loadingData ||
                              loadingParams ||
                              loadingCalculate ||
                              resultsTable.data.length
                            }
                            type="checkbox"
                            value={exposureVar.log2}
                            checked={exposureVar.log2}
                            onChange={() =>
                              mergeState({
                                exposureVar: { log2: !exposureVar.log2 },
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
                          !exposureVar.collapseOptions.length
                        }
                        id="collapse2"
                        label="Collapse"
                        value={
                          exposureVar.collapseOptions.length
                            ? exposureVar.collapse
                            : 'None'
                        }
                        options={exposureVar.collapseOptions}
                        onChange={(e) =>
                          mergeState({ exposureVar: { collapse: e } })
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
                    disabled={
                      loadingData ||
                      loadingParams ||
                      loadingCalculate ||
                      resultsTable.data.length
                    }
                    id="testType"
                    label=""
                    value={testType}
                    options={['lm', 'glm']}
                    onChange={(e) => mergeState({ testType: e })}
                  />
                </fieldset>
              </Col>
              <Col md="auto" className="d-flex">
                <Button
                  disabled={
                    loadingData ||
                    loadingParams ||
                    loadingCalculate ||
                    resultsTable.data.length
                  }
                  className="w-100 align-self-center"
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
                  ]
                    .reduce(reducer, [])
                    .map((col) =>
                      typeof resultsTable.data[0][col.Header] == 'number'
                        ? { ...col, sortType: 'basic' }
                        : col
                    )}
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
              </Row>
              <LoadingOverlay active={loadingRecalculate} />
              {plotPath && (
                <Plot
                  className="p-3 border rounded"
                  downloadName={plotPath.split('/').slice(-1)[0]}
                  plotPath={`api/results/${plotPath}`}
                  txtPath={`api/results/${dataPath}`}
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
