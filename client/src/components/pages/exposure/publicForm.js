import React, { useEffect, useState } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import Select from '../../controls/select/select';
import { useSelector, useDispatch } from 'react-redux';
import { actions as exposureActions } from '../../../services/store/exposure';
import { actions as modalActions } from '../../../services/store/modal';
import { getJSON } from '../../../services/utils';

const actions = { ...exposureActions, ...modalActions };
const { Group, Check } = Form;

export default function PublicForm({
  calculate,
  handleReset,
  handleStudy,
  handleStrategy,
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
    exposureSignature,
    studyOptions,
    strategyOptions,
    cancerOptions,
    rsSetOptions,
    study,
    strategy,
    cancer,
    rsSet,
    gettingSignatureNames,
    gettingSampleNames,
    useCancerType,
  } = useSelector((state) => state.exposure.exposureState);

  const [queryNames, setQuery] = useState(false);

  // populate controls on inital render
  useEffect(() => {
    if (!exposureSignature.length) populateControls();
  }, []);

  // call calculate after receiving signature and sample name options
  useEffect(() => {
    if (!submitted && queryNames) calculate();
  }, [queryNames, submitted]);

  async function populateControls() {
    try {
      let [
        exposureCancer,
        exposureSignature,
        signatureNames,
      ] = await Promise.all([
        getJSON('Others/json/Exploring-Exposure-cancertype.json'),
        getJSON('Others/json/Exploring-Exposure.json'),
        getJSON('Others/json/Signature_name.json'),
      ]);

      populateExposureExp(exposureCancer, exposureSignature, signatureNames);
    } catch (err) {
      mergeError(err.message);
    }
  }

  async function populateExposureExp(
    exposureCancer,
    exposureSignature,
    signatureNames
  ) {
    mergeState({ loading: true });

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
        exposureCancer
          .filter((data) => data.Study == study && data.Dataset == strategy)
          .map((data) => data.Cancer_Type)
      ),
    ];
    const cancer = 'Lung-AdenoCA'; // default

    const signatureNameOptions = [
      ...new Set(
        signatureNames
          .filter((row) => row.Signature_set_name == rsSet)
          .map((row) => row.Signature_name)
      ),
    ];

    mergeState({
      study: study,
      studyOptions: studyOptions,
      strategy: strategy,
      strategyOptions: strategyOptions,
      cancer: cancer,
      cancerOptions: cancerOptions,
      rsSet: rsSet,
      rsSetOptions: rsSetOptions,
      signatureNameOptions: signatureNameOptions,
      loading: false,
      exposureCancer,
      exposureSignature,
      signatureNames,
    });
  }

  // get sample name options filtered by cancer type
  async function getSampleNames() {
    try {
      const { stdout, output } = await (
        await fetch(`api/explorationData`, {
          method: 'POST',
          headers: {
            Accept: 'application/json',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            fn: 'getSampleNames',
            args: {
              study: study,
              strategy: strategy,
              rsSet: rsSet,
              cancerType: cancer,
            },
          }),
        })
      ).json();

      if (output.data.length) {
        mergeState({
          publicSampleOptions: output.data,
        });
        // mergeMsIndividual({ sample: output.data[0] });
      } else console.log('No Sample Names Found');
    } catch (err) {
      mergeError(err.message);
    }
  }

  // get signature name options filtered by cancer type
  async function getSignatureNames() {
    if (useCancerType) {
      try {
        const { stdout, output } = await (
          await fetch(`api/explorationData`, {
            method: 'POST',
            headers: {
              Accept: 'application/json',
              'Content-Type': 'application/json',
            },
            body: JSON.stringify({
              fn: 'getSignatureNames',
              args: {
                study: study,
                strategy: strategy,
                rsSet: rsSet,
                cancerType: cancer,
              },
            }),
          })
        ).json();

        if (output.data.length)
          mergeState({
            signatureNameOptions: output.data,
          });
        else console.log('No Signature Names Found');
      } catch (err) {
        mergeError(err.message);
      }
    }
  }

  //   get signature and sample names. useEffect will call main calculate function
  async function handleCalculate() {
    mergeState({ loading: true });
    await Promise.all([getSampleNames(), getSignatureNames()]);
    mergeState({ loading: false });
    setQuery(true);
  }

  return (
    <Form>
      <Row>
        <Col>
          <Group>
            <Select
              disabled={
                loading ||
                submitted ||
                gettingSignatureNames ||
                gettingSampleNames
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
                loading ||
                submitted ||
                gettingSignatureNames ||
                gettingSampleNames
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
                loading ||
                submitted ||
                gettingSignatureNames ||
                gettingSampleNames
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
                loading ||
                submitted ||
                gettingSignatureNames ||
                gettingSampleNames
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
        <Col>
          <Group controlId="toggleCancerType">
            <Check
              disabled={
                loading ||
                submitted ||
                gettingSampleNames ||
                gettingSignatureNames
              }
              type="checkbox"
              label="Cancer Type Only"
              value={useCancerType}
              checked={useCancerType}
              onChange={(e) => {
                if (!useCancerType == false) handleSet(rsSet);

                mergeState({
                  useCancerType: !useCancerType,
                });
              }}
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
              setQuery(false);
              handleReset();
            }}
          >
            Reset
          </Button>
        </Col>
        <Col lg="6">
          <Button
            disabled={loading || submitted}
            className="w-100"
            variant="primary"
            onClick={() => handleCalculate()}
          >
            Calculate
          </Button>
        </Col>
      </Row>
    </Form>
  );
}
