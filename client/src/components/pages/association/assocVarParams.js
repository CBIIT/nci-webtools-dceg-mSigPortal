import React, { useEffect } from 'react';
import {
  Form,
  Row,
  Col,
  Button,
  OverlayTrigger,
  Popover,
} from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faInfoCircle, faMinus } from '@fortawesome/free-solid-svg-icons';
import Select from '../../controls/select/select';

const { Group, Label, Check, Control } = Form;

export default function AssocVarParams({
  hostState,
  paramState,
  mergeState,
  handleLoadParameters,
  remove = false,
  last = true,
}) {
  const { loadingData, assocVarData } = useSelector(
    (state) => state.association.associationState
  );

  const { loadingParams, loadingCalculate } = hostState;

  const {
    source = '',
    type = '',
    tmpName = '',
    sourceOptions = [],
    typeOptions = [],
    nameOptions = [],
    filter = '',
    log2 = false,
    collapse = '',
    collapseOptions = [],
  } = paramState;

  // populate controls
  useEffect(() => {
    if (assocVarData.length && !source) {
      const sourceOptions = [
        ...new Set(assocVarData.map((row) => row.data_source)),
      ];
      const source = sourceOptions[0];
      console.log(source, sourceOptions);

      mergeState({ source, sourceOptions });
    }
  }, [assocVarData]);

  // filter dropdown options on change
  useEffect(() => {
    if (source) handleSource();
  }, [source]);
  useEffect(() => {
    if (type) handleType();
  }, [type]);

  function handleSource() {
    const typeOptions = [
      ...new Set(
        assocVarData
          .filter((row) => row.data_source == source)
          .map((row) => row.data_type)
      ),
    ];

    mergeState({ type: typeOptions[0], typeOptions });
  }

  function handleType() {
    const nameOptions = [
      ...new Set(
        assocVarData
          .filter((row) => row.data_type == type)
          .map((row) => row.variable_name)
      ),
    ];

    mergeState({ tmpName: nameOptions[0], nameOptions });
  }

  const popoverInfo = (text) => (
    <Popover>
      <Popover.Content>{text}</Popover.Content>
    </Popover>
  );

  return (
    <div>
      <Row className="justify-content-center mt-3">
        <Col md="auto" className="d-flex">
          {remove ? (
            <Button
              className="text-danger mb-3"
              variant="link"
              onClick={remove}
              title="Add Plot"
              style={{ textDecoration: 'none' }}
            >
              <span className="text-nowrap" title="Remove Variable">
                <FontAwesomeIcon icon={faMinus} /> Remove
              </span>
            </Button>
          ) : (
            <span style={{ width: '101.867px' }} />
          )}
        </Col>
        <Col md="auto">
          <Select
            disabled={loadingData || loadingParams || loadingCalculate}
            id="source"
            label="Variable Source"
            value={source}
            options={sourceOptions}
            onChange={(e) => mergeState({ source: e })}
          />
        </Col>
        <Col md="auto" lg="3" xl="2">
          <Select
            disabled={loadingData || loadingParams || loadingCalculate}
            id="type"
            label="Data Type"
            value={type}
            options={typeOptions}
            onChange={(e) => mergeState({ type: e })}
          />
        </Col>
        <Col md="auto" lg="3" xl="2">
          <Select
            disabled={loadingData || loadingParams || loadingCalculate}
            id="assocVariable"
            label="Variant Name"
            value={tmpName}
            options={nameOptions}
            onChange={(e) => mergeState({ tmpName: e })}
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
                  <Group controlId="threshold" className="d-flex mb-0">
                    <Label className="mr-2 font-weight-normal">Threshold</Label>
                    <Control
                      disabled={
                        loadingData || loadingParams || loadingCalculate
                      }
                      value={filter}
                      placeholder={'Optional'}
                      style={{ width: '90px' }}
                      onChange={(e) => mergeState({ filter: e.target.value })}
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
                      value={log2}
                      checked={log2}
                      onChange={() => mergeState({ log2: !log2 })}
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
                </OverlayTrigger>
                <Select
                  className="d-inline-flex mb-0"
                  disabled={
                    loadingData ||
                    loadingParams ||
                    loadingCalculate ||
                    !collapseOptions.length
                  }
                  id="collapse1"
                  label="Collapse"
                  labelClass="mr-2 font-weight-normal"
                  value={collapseOptions.length ? collapse : 'None'}
                  options={collapseOptions}
                  onChange={(e) => mergeState({ collapse: e })}
                />
              </Col>
            </Row>
          </fieldset>
        </Col>
        <Col md="auto" className="d-flex">
          {last ? (
            <Button
              disabled={
                loadingData || loadingParams || loadingCalculate || !source
              }
              className="w-100 align-self-center"
              variant="primary"
              onClick={() => handleLoadParameters()}
            >
              Load Data
            </Button>
          ) : (
            <span style={{ width: '120px' }} />
          )}
        </Col>
      </Row>
    </div>
  );
}