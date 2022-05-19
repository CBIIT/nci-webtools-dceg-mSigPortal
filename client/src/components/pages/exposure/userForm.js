import React, { useState, useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faFolderMinus } from '@fortawesome/free-solid-svg-icons';
import CustomSelect from '../../controls/select/select';
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
    projectID,
    userNameOptions,
    userSampleOptions,
  } = useSelector((state) => state.exposure.exposureState);

  const [exposureFileObj, setExposure] = useState(new File([], ''));
  const [matrixFileObj, setMatrix] = useState(new File([], ''));
  const [signatureFileObj, setSignature] = useState(new File([], ''));

  const [checkValid, setCheckValid] = useState(false);

  // call calculate after receiving signature and sample name options
  useEffect(() => {
    if (
      !submitted &&
      !loading &&
      userNameOptions.length &&
      userSampleOptions.length
    )
      calculate('all', projectID);
  }, [userNameOptions, userSampleOptions, loading]);

  function validateFiles() {
    setCheckValid(true);
    return usePublicSignature
      ? exposureFileObj.size && matrixFileObj.size
      : exposureFileObj.size && matrixFileObj.size && signatureFileObj.size;
  }

  async function loadExample(type) {
    const filepath = `assets/exampleInput/Sherlock_SBS96_${type}.txt`;
    const filename = filepath.split('/').slice(-1)[0];
    if (`${type}File` != filename) {
      if (type == 'exposure') {
        setExposure(new File([await (await fetch(filepath)).blob()], filename));
        mergeState({ [`${type}File`]: filename });
      } else if (type == 'matrix') {
        setMatrix(new File([await (await fetch(filepath)).blob()], filename));
        mergeState({ [`${type}File`]: filename });
      } else if (type == 'signature') {
        setSignature(
          new File([await (await fetch(filepath)).blob()], filename)
        );
        mergeState({ [`${type}File`]: filename });
      }
    }
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

  //  get signature and sample names. useEffect will call main calculate function
  async function handleCalculate() {
    mergeState({ loading: true });
    try {
      const { projectID, exposureData } = await handleUpload();
      // get signature name options, ignore sample key
      const nameOptions = Object.keys(exposureData[0]).filter(
        (key) => key != 'Samples'
      );
      const sampleOptions = [
        ...new Set(exposureData.map(({ Samples }) => Samples)),
      ];

      dispatch(
        actions.mergeExposure({
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
        })
      );
    } catch (err) {
      mergeError(err);
    }
    mergeState({ loading: false });
  }

  return (
    <Form>
      <Row>
        <Col>
          <Group>
            <Label>Upload Exposure File</Label>
            <Row className="m-0">
              <Col lg="6" className="p-0">
                <Button
                  className="p-0 font-14"
                  disabled={submitted}
                  variant="link"
                  href={'assets/exampleInput/Sherlock_SBS96_exposure.txt'}
                  download
                >
                  Download Example
                </Button>
              </Col>
              <Col lg="6" className="p-0 d-flex">
                <Button
                  className={`p-0 ml-auto font-14`}
                  disabled={submitted || loading.active}
                  variant="link"
                  type="button"
                  onClick={() => loadExample('exposure')}
                >
                  Load Example
                </Button>
              </Col>
            </Row>
            <div className="d-flex">
              <Form.File
                disabled={loading || submitted}
                id="uploadExposure"
                label={exposureFileObj.name || 'Exposure File'}
                title={exposureFileObj.name || 'Upload Exposure File'}
                value={''}
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
                  disabled={loading || submitted}
                  onClick={() => {
                    setExposure(new File([], ''));
                    mergeState({ exposureFile: '' });
                  }}
                >
                  <FontAwesomeIcon icon={faFolderMinus} size="lg" />
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
            <Row className="m-0">
              <Col lg="6" className="p-0">
                <Button
                  className="p-0 font-14"
                  disabled={submitted}
                  variant="link"
                  href={'assets/exampleInput/Sherlock_SBS96_matrix.txt'}
                  download
                >
                  Download Example
                </Button>
              </Col>
              <Col lg="6" className="p-0 d-flex">
                <Button
                  className={`p-0 ml-auto font-14`}
                  disabled={submitted || loading.active}
                  variant="link"
                  type="button"
                  onClick={() => loadExample('matrix')}
                >
                  Load Example
                </Button>
              </Col>
            </Row>
            <div className="d-flex">
              <Form.File
                disabled={loading || submitted}
                id="uploadMatrix"
                label={matrixFileObj.name || 'Matrix File'}
                title={matrixFileObj.name || 'Upload Matrix File'}
                value={''}
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
                  disabled={loading || submitted}
                  onClick={() => {
                    setMatrix(new File([], ''));
                    mergeState({ matrixFile: '' });
                  }}
                >
                  <FontAwesomeIcon icon={faFolderMinus} size="lg" />
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
                <CustomSelect
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
                <CustomSelect
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
              <Row className="m-0">
                <Col lg="6" className="p-0">
                  <Button
                    className="p-0 font-14"
                    disabled={submitted}
                    variant="link"
                    href={'assets/exampleInput/Sherlock_SBS96_signature.txt'}
                    download
                  >
                    Download Example
                  </Button>
                </Col>
                <Col lg="6" className="p-0 d-flex">
                  <Button
                    className={`p-0 ml-auto font-14`}
                    disabled={submitted || loading.active}
                    variant="link"
                    type="button"
                    onClick={() => loadExample('signature')}
                  >
                    Load Example
                  </Button>
                </Col>
              </Row>
              <div className="d-flex">
                <Form.File
                  disabled={loading || submitted}
                  id="uploadSignature"
                  label={signatureFileObj.name || 'Signature File'}
                  title={signatureFileObj.name || 'Upload Signature File'}
                  value={''}
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
                    disabled={loading || submitted}
                    onClick={() => {
                      setSignature(new File([], ''));
                      mergeState({ signatureFile: '' });
                    }}
                  >
                    <FontAwesomeIcon icon={faFolderMinus} size="lg" />
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
            <CustomSelect
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
