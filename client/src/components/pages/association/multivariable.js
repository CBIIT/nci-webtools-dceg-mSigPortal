import React, { useEffect, useRef } from 'react';
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

  async function handleLoadParameters() {
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
            fn: 'loadCollapse',
            args: {
              study,
              strategy,
              rsSet,
              cancer,
              source: associationVars.source,
              type: associationVars.type,
              assocName: associationVars.tmpName,
              expName: exposureVar.name,
            },
          }),
        })
      ).json();

      const { collapseVar1, collapseVar2 } = collapseData;

      mergeState({
        associationVars: {
          name: associationVars.tmpName,
          collapseOptions: collapseVar1 || [],
        },
        exposureVar: {
          name: expVarList[0],
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
              associationVars: (() => {
                const {
                  sourceOptions,
                  typeOptions,
                  nameOptions,
                  tmpName,
                  collapseOptions,
                  ...params
                } = associationVars;
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
    let newParams = associationVars.slice();
    newParams.push({});
    mergeState({ associationVars: newParams });
  }
  function removeParam(index) {
    let newParams = associationVars.slice();
    newParams.splice(index, 1);
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
        <LoadingOverlay active={loadingParams} />
        <p className="text-center">
          Select the following variables for analysis
        </p>

        {associationVars.map((paramState, index) => (
          <AssocVarParams
            hostState={multivariable}
            paramState={paramState}
            mergeState={(e) => {
              let newParams = associationVars.slice();
              newParams[index] = { ...newParams[index], ...e };
              mergeState({ associationVars: newParams });
            }}
            handleLoadParameters={() => {}} //handleLoadParameters}
            remove={index != 0 ? () => removeParam(index) : false}
            last={index == associationVars.length - 1}
          />
        ))}
        <Row className="mt-3">
          <Col md="auto" className="d-flex">
            <Button
              className="ml-auto"
              variant="link"
              onClick={() => addParam()}
              title="Add Plot"
              style={{ textDecoration: 'none' }}
            >
              <span className="text-nowrap" title="Add Assocation Variable">
                <FontAwesomeIcon icon={faPlus} /> Add Association Variable
              </span>
            </Button>
          </Col>
        </Row>
      </div>
      <div className="mb-3">
        <h5 className="separator">Parameters</h5>

        <LoadingOverlay active={loadingCalculate} />
        {associationVars.name && exposureVar.name ? (
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
                              loadingData || loadingParams || loadingCalculate
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
                              loadingData || loadingParams || loadingCalculate
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
                    disabled={loadingData || loadingParams || loadingCalculate}
                    id="testType"
                    label=""
                    value={testType}
                    options={['nonparametric', 'parametric']}
                    onChange={(e) => mergeState({ testType: e })}
                  />
                </fieldset>
              </Col>
              <Col md="auto" className="d-flex">
                <Button
                  disabled={loadingData || loadingParams || loadingCalculate}
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
                <Col md="auto">
                  <Group controlId="xlab">
                    <Label>Variable Title</Label>
                    <Control
                      value={xlab}
                      placeholder={associationVars.name}
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
