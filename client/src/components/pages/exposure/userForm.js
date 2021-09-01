import React, { useState } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faFolderMinus } from '@fortawesome/free-solid-svg-icons';
import Select from '../../controls/select/select';
import { useSelector, useDispatch } from 'react-redux';
import { actions as exposureActions } from '../../../services/store/exposure';
import { actions as modalActions } from '../../../services/store/modal';

const actions = { ...exposureActions, ...modalActions };
const { Group, Check, Label } = Form;

export default function PublicForm({
  calculate,
  handleReset,
  handleStudy,
  handleSet,
}) {
  const dispatch = useDispatch();
  const mergeState = async (state) =>
    await dispatch(actions.mergeExposure({ exposureState: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const {
    submitted,
    loading,
    studyOptions,
    rsSetOptions,
    study,
    rsSet,
    usePublicSignature,
    genome,
    genomeOptions,
  } = useSelector((state) => state.exposure.exposureState);

  const [exposureFileObj, setExposure] = useState(new File([], ''));
  const [matrixFileObj, setMatrix] = useState(new File([], ''));
  const [signatureFileObj, setSignature] = useState(new File([], ''));

  const [checkValid, setCheckValid] = useState(false);

  function validateFiles() {
    setCheckValid(true);
    return usePublicSignature
      ? exposureFileObj.size && matrixFileObj.size
      : exposureFileObj.size && matrixFileObj.size && signatureFileObj.size;
  }

  async function handleUpload() {
    return new Promise(async (resolve, reject) => {
      if (validateFiles()) {
        try {
          const data = new FormData();
          data.append('exposureFile', exposureFileObj);
          data.append('matrixFile', matrixFileObj);
          if (!usePublicSignature)
            data.append('signatureFile', signatureFileObj);
          // if (variableFileObj.size)
          //   dasel'variableFile', variableFileObj);
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
            const { projectID, exposurePath } = await response.json();

            const exposureData = await (
              await fetch('api/getSignaturesUser', {
                method: 'POST',
                headers: {
                  Accept: 'application/json',
                  'Content-Type': 'application/json',
                },
                body: JSON.stringify({
                  path: exposurePath,
                }),
              })
            ).json();
            await mergeState({ projectID });
            resolve({ projectID, exposureData });
          }
        } catch (err) {
          mergeError(err.message);
          reject(err);
        }
      } else {
        reject('Missing required files');
      }
    });
  }

  async function handleCalculate() {
    try {
      const { projectID, exposureData } = await handleUpload();
      // get signature name options, ignore sample key
      const nameOptions = exposureData.columns.filter(
        (key) => key != 'Samples'
      );

      const sampleOptions = [
        ...new Set(exposureData.map(({ Samples }) => Samples)),
      ];

      const params = {
        exposureState: {
          projectID: projectID,
          userNameOptions: nameOptions,
          userSampleOptions: sampleOptions,
        },
        msIndividual: { sample: sampleOptions[0] },
        msAssociation: {
          signatureName1: nameOptions[0],
          signatureName2: nameOptions[1],
        },
        msBurden: { signatureName: nameOptions[0] },
      };

      await dispatch(actions.mergeExposure(params));
      await calculate('all', projectID, params);
    } catch (err) {
      mergeError(err);
    }
  }

  return (
    <Form>
      <Row>
        <Col>
          <Group>
            <Label>Upload Exposure File</Label>
            <div className="d-flex">
              <Form.File
                disabled={loading || submitted}
                id="uploadExposure"
                label={exposureFileObj.name || 'Exposure File'}
                accept=".txt"
                isInvalid={checkValid && !exposureFileObj.size}
                feedback="Upload an exposure file"
                onChange={(e) => {
                  if (e.target.files.length) {
                    setExposure(e.target.files[0]);
                    mergeState({ exposureFile: e.target.files[0].name });
                  }
                }}
                custom
              />
              {exposureFileObj.size > 0 && (
                <Button
                  className="ml-1"
                  size="sm"
                  title="Remove"
                  variant="danger"
                  onClick={() => {
                    setExposure(new File([], ''));
                    mergeState({ exposureFile: '' });
                  }}
                >
                  <FontAwesomeIcon icon={faFolderMinus} size="2x" />
                </Button>
              )}
            </div>
          </Group>
        </Col>
      </Row>
      <Row>
        <Col>
          <Group>
            <Label>Upload Matrix File</Label>
            <div className="d-flex">
              <Form.File
                disabled={loading || submitted}
                id="uploadMatrix"
                label={matrixFileObj.name || 'Matrix File'}
                accept=".txt"
                isInvalid={checkValid && !matrixFileObj.size}
                feedback="Upload a matrix file"
                onChange={(e) => {
                  if (e.target.files.length) {
                    setMatrix(e.target.files[0]);
                    mergeState({ matrixFile: e.target.files[0].name });
                  }
                }}
                custom
              />
              {matrixFileObj.size > 0 && (
                <Button
                  className="ml-1"
                  size="sm"
                  title="Remove"
                  variant="danger"
                  onClick={() => {
                    setMatrix(new File([], ''));
                    mergeState({ matrixFile: '' });
                  }}
                >
                  <FontAwesomeIcon icon={faFolderMinus} size="2x" />
                </Button>
              )}
            </div>
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
                  disabled={loading || submitted || !studyOptions.length}
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
                  disabled={loading || submitted || !rsSetOptions.length}
                  id="exposureSignatureSet"
                  label="Reference Signature Set"
                  value={rsSet}
                  options={rsSetOptions}
                  onChange={handleSet}
                />
              </Group>
            </Col>
          </Row>
        </div>
      ) : (
        <Row>
          <Col>
            <Group>
              <Label>Upload Signature Data</Label>
              <div className="d-flex">
                <Form.File
                  disabled={loading || submitted}
                  id="uploadSignature"
                  label={signatureFileObj.name || 'Signature File'}
                  accept=".txt"
                  isInvalid={checkValid && !signatureFileObj.size}
                  feedback="Upload a signature file"
                  onChange={(e) => {
                    if (e.target.files.length) {
                      setSignature(e.target.files[0]);
                      mergeState({ signatureFile: e.target.files[0].name });
                    }
                  }}
                  custom
                />
                {signatureFileObj.size > 0 && (
                  <Button
                    className="ml-1"
                    size="sm"
                    title="Remove"
                    variant="danger"
                    onClick={() => {
                      setSignature(new File([], ''));
                      mergeState({ signatureFile: '' });
                    }}
                  >
                    <FontAwesomeIcon icon={faFolderMinus} size="2x" />
                  </Button>
                )}
              </div>
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
        <Col lg="6">
          <Button
            disabled={loading}
            className="w-100 mb-3"
            variant="secondary"
            onClick={() => {
              setCheckValid(false);
              setExposure(new File([], ''));
              setMatrix(new File([], ''));
              setSignature(new File([], ''));
              handleReset();
            }}
          >
            Reset
          </Button>
        </Col>
        <Col lg="6">
          <Button
            disabled={loading}
            className="w-100"
            variant="primary"
            onClick={() => {
              if (validateFiles()) handleCalculate();
            }}
          >
            Calculate
          </Button>
        </Col>
      </Row>
    </Form>
  );
}
