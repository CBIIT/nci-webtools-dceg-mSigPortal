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
import { faInfoCircle } from '@fortawesome/free-solid-svg-icons';
import { actions as associationActions } from '../../../services/store/association';
import { actions as modalActions } from '../../../services/store/modal';
import CustomSelect from '../../controls/select/select-old';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import AssocVarParams from './assocVarParams';
import SvgContainer from '../../controls/svgContainer/svgContainer';
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
  } = useSelector((state) => state.association.main);

  const {
    loadingParams,
    loadingCalculate,
    loadingRecalculate,
    error,
    loadError,
    id,
    plotPath,
    dataPath,
    assocTablePath,
    signature,
    signatureOptions,
    testType,
    xlab,
    ylab,
    associationVar,
    exposureVar,
    resultsTable,
  } = useSelector((state) => state.association.univariable);

  const [invalidAssocFilter, setInvalidAssocFilter] = useState(false);
  const [invalidExpFilter, setInvalidExpFilter] = useState(false);

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
  async function handleLoadData() {
    if (associationVar.filter && isNaN(associationVar.filter)) {
      setInvalidAssocFilter(true);
    } else {
      setInvalidAssocFilter(false);

      mergeState({ loadingParams: true, error: false });
      try {
        const {
          sessionId,
          stdout,
          output: collapseData,
        } = await (
          await fetch(`web/associationWrapper`, {
            method: 'POST',
            headers: {
              Accept: 'application/json',
              'Content-Type': 'application/json',
            },
            body: JSON.stringify({
              fn: 'loadCollapse',
              args: {
                study: study.value,
                strategy: strategy.value,
                rsSet: rsSet.value,
                cancer: cancer.value,
                source: associationVar.source,
                type: associationVar.type,
                assocName: associationVar.tmpName,
                expName: exposureVar.name,
              },
            }),
          })
        ).json();
        const { collapseVar1, collapseVar2, error, uncaughtError } =
          collapseData;

        if (error || uncaughtError) {
          mergeState({
            loadError:
              error ||
              'An error has occured. Please review your input and try again. If the issue persists, please contact us: NCImSigPortalWebAdmin@mail.nih.gov',
          });
        } else {
          mergeState({
            loadError: '',
            associationVar: {
              name: associationVar.tmpName,
              collapseOptions: Array.isArray(collapseVar1) ? collapseVar1 : [],
              id: sessionId,
            },
            // exposureVar: {
            //   name: expVarList[0],
            //   collapseOptions: Array.isArray(collapseVar2) ? collapseVar2 : [],
            // },
          });
        }
      } catch (error) {
        mergeError(error);
      }
      mergeState({ loadingParams: false });
    }
  }
  function handleReset() {
    setInvalidAssocFilter(false);
    setInvalidExpFilter(false);
    mergeState({
      associationVar: {
        name: '',
        filter: '',
        log2: false,
        collapse: '',
        collapseOptions: [],
      },
      id: '',
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
    if (exposureVar.filter && isNaN(exposureVar.filter)) {
      setInvalidExpFilter(true);
    } else {
      setInvalidExpFilter(false);

      mergeState({
        loadingCalculate: signature ? false : true,
        loadingRecalculate: signature ? true : false,
        error: false,
        plotPath: '',
        dataPath: '',
      });
      try {
        const { sessionId, stdout, output } = await (
          await fetch(`web/associationWrapper`, {
            method: 'POST',
            headers: {
              Accept: 'application/json',
              'Content-Type': 'application/json',
            },
            body: JSON.stringify({
              fn: 'univariable',
              id,
              args: {
                study: study.value,
                strategy: strategy.value,
                rsSet: rsSet.value,
                cancer: cancer.value,
                testType: testType,
                signature: signature,
                xlab: xlab || associationVar.name,
                ylab: ylab || exposureVar.name,
                associationVar: (() => {
                  const {
                    sourceOptions,
                    typeOptions,
                    nameOptions,
                    tmpName,
                    collapseOptions,
                    ...params
                  } = associationVar;
                  return params;
                })(),
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
          mergeState({
            error:
              error ||
              uncaughtError ||
              'An error has occured. Please review your input and try again. If the issue persists, please contact us: NCImSigPortalWebAdmin@mail.nih.gov',
          });
        } else {
          mergeState({
            id,
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
          downloadLink={assocFullDataPath}
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

          <AssocVarParams
            hostState={useSelector((state) => state.association.univariable)}
            paramState={associationVar}
            mergeState={(e) => mergeState({ associationVar: e })}
            invalidFilter={invalidAssocFilter}
          />
          <Row
            className="mx-auto mt-3 justify-content-between"
            style={{ maxWidth: '1720px' }}
          >
            <Col md="auto" lg="auto">
              <CustomSelect
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
            <Col md="auto">
              <Button
                disabled={loadingData || loadingParams || loadingCalculate}
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
                  !associationVar.source ||
                  associationVar.name
                }
                variant="primary"
                onClick={() => handleLoadData()}
              >
                Load Data
              </Button>
            </Col>
          </Row>
          {loadError && (
            <p className="p-3 d-flex justify-content-center text-danger">
              {loadError}
            </p>
          )}
        </div>
      </div>
      <div className="mb-3">
        <h5 className="separator">Parameters</h5>

        <LoadingOverlay active={loadingCalculate} />
        {associationVar.name && exposureVar.name ? (
          <>
            <p className="text-center">
              Select the following filtering and method for analysis
            </p>
            <Row className="justify-content-center">
              <Col lg="auto">
                <fieldset className="border rounded p-2">
                  <legend className="font-weight-bold">
                    Signature Exposure Filtering
                  </legend>
                  <Row>
                    <Col md="auto">
                      <Group controlId="exposureVar-threshold" className="mb-0">
                        <div className="d-flex">
                          <Label className="mr-2 font-weight-normal">
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
                                  style={{ verticalAlign: 'baseline' }}
                                />
                              </Button>
                            </OverlayTrigger>
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
                            isInvalid={invalidExpFilter}
                          />
                        </div>
                        {invalidExpFilter && (
                          <Form.Control.Feedback
                            className="d-block"
                            type="invalid"
                          >
                            Enter a numeric threshold value
                          </Form.Control.Feedback>
                        )}
                      </Group>
                    </Col>
                    <Col md="auto">
                      <Group controlId="log2-2" className="d-flex mb-0">
                        <Label className="mr-2 font-weight-normal">
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
                      <CustomSelect
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
                  <CustomSelect
                    aria-label="Method"
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
                    options={[
                      'nonparametric',
                      'parametric',
                      'robust',
                      'bayes',
                      'skit',
                    ]}
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
                  title="Variable Association"
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
                  downloadLink={id + assocTablePath}
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
                  <CustomSelect
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
                      placeholder={associationVar.name}
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
                      placeholder={exposureVar.name}
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
                <SvgContainer
                  className="p-3 border rounded"
                  downloadName={plotPath.split('/').slice(-1)[0]}
                  plotPath={`web/data/${plotPath}`}
                  txtPath={`web/data/${dataPath}`}
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
