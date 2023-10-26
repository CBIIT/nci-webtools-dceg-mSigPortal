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
import {
  faInfoCircle,
  faMinus,
  faPlus,
} from '@fortawesome/free-solid-svg-icons';
import CustomSelect from '../../controls/select/select-old';

const { Group, Label, Check, Control } = Form;

export default function AssocVarParams({
  index = 0,
  hostState,
  paramState,
  mergeState,
  remove = false,
  add = false,
  warnLimit = false,

  duplicates = [],
  invalidFilter,
}) {
  const { loadingData, assocVarData } = useSelector(
    (state) => state.association.main
  );

  const { loadingParams, loadingCalculate, associationVars } = hostState;

  const {
    name,
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

  const getTypeOptions = (source) => [
    ...new Set(
      assocVarData
        .filter((row) => row.data_source == source)
        .map((row) => row.data_type)
    ),
  ];

  const getNameOptions = (source, type) => [
    ...new Set(
      assocVarData
        .filter((row) => row.data_source == source && row.data_type == type)
        .map((row) => row.variable_name)
    ),
  ];

  // populate controls
  useEffect(() => {
    if (assocVarData.length && !source) {
      const sourceOptions = [
        ...new Set(assocVarData.map((row) => row.data_source)),
      ];
      const source = sourceOptions[0];
      const typeOptions = getTypeOptions(source);
      const type = typeOptions[0];
      const nameOptions = getNameOptions(source, type);

      mergeState({
        source,
        sourceOptions,
        type,
        typeOptions,
        nameOptions,
        tmpName: nameOptions[0],
      });
    }
  }, [source]);

  function handleSource(source) {
    const typeOptions = getTypeOptions(source);
    const nameOptions = getNameOptions(source, typeOptions[0]);

    mergeState({
      source,
      type: typeOptions[0],
      typeOptions,
      tmpName: nameOptions[0],
      nameOptions,
    });
  }

  function handleType(type) {
    const nameOptions = getNameOptions(source, type);

    mergeState({
      type,
      tmpName: nameOptions[0],
      nameOptions,
    });
  }

  const popoverInfo = (text) => (
    <Popover>
      <Popover.Content>{text}</Popover.Content>
    </Popover>
  );

  return (
    <div>
      <Row
        className={`justify-content-center mt-3 border rounded p-2 ${
          duplicates.indexOf(index) > -1 ? 'border-danger' : ''
        }`}
      >
        <Col md="auto">
          <CustomSelect
            disabled={loadingData || loadingParams || loadingCalculate || name}
            id={'source-' + index}
            label="Variable Source"
            value={source}
            options={sourceOptions}
            onChange={(e) => handleSource(e)}
          />
        </Col>
        <Col md="auto">
          <CustomSelect
            disabled={loadingData || loadingParams || loadingCalculate || name}
            id={'type-' + index}
            label="Data Type"
            value={type}
            options={typeOptions}
            onChange={(e) => handleType(e)}
          />
        </Col>
        <Col md="auto">
          <CustomSelect
            disabled={loadingData || loadingParams || loadingCalculate || name}
            id={'assocVariable-' + index}
            label="Variable Name"
            value={tmpName}
            options={nameOptions}
            onChange={(e) => mergeState({ tmpName: e })}
          />
        </Col>

        <Col md="auto">
          <fieldset className="border rounded p-2">
            <legend className="font-weight-bold">Variable Filtering</legend>
            <Row>
              <Col md="auto">
                <Group controlId={'threshold-' + index} className="mb-0">
                  <div className="d-flex">
                    <Label className="mr-2 font-weight-normal">
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
                            style={{ verticalAlign: 'baseline' }}
                          />
                        </Button>
                      </OverlayTrigger>
                      Threshold
                    </Label>
                    <Control
                      disabled={
                        loadingData || loadingParams || loadingCalculate || name
                      }
                      value={filter}
                      placeholder={'Optional'}
                      style={{ width: '90px' }}
                      onChange={(e) => mergeState({ filter: e.target.value })}
                      isInvalid={invalidFilter}
                    />
                  </div>
                  {invalidFilter && (
                    <Form.Control.Feedback className="d-block" type="invalid">
                      Enter a numeric threshold value
                    </Form.Control.Feedback>
                  )}
                </Group>
              </Col>
              <Col md="auto">
                <Group controlId={'log2-1-' + index} className="d-flex mb-0">
                  <Label className="mr-2 font-weight-normal">
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
                    </OverlayTrigger>
                    log<sub>2</sub>
                  </Label>
                  <Check>
                    <Check.Input
                      disabled={
                        loadingData || loadingParams || loadingCalculate || name
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
                </OverlayTrigger>{' '}
                <CustomSelect
                  className="d-inline-flex mb-0"
                  disabled={
                    loadingData ||
                    loadingParams ||
                    loadingCalculate ||
                    !collapseOptions.length ||
                    name
                  }
                  id={'collapse1-' + index}
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
          {remove ? (
            <Button
              disabled={name}
              className="text-danger mr-auto"
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
          {add ? (
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
                disabled={
                  loadingData ||
                  loadingParams ||
                  loadingCalculate ||
                  associationVars[0].name
                }
                variant="link"
                onClick={add}
                title="Add Plot"
                style={{ textDecoration: 'none' }}
              >
                <span className="text-nowrap" title="Add Association Variable">
                  <FontAwesomeIcon icon={faPlus} /> Add Association Variable
                </span>
              </Button>
            </OverlayTrigger>
          ) : (
            <span style={{ width: '101.867px' }} />
          )}
        </Col>
      </Row>
    </div>
  );
}
