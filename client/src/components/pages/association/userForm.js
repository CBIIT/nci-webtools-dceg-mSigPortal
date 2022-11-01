import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import CustomSelect from '../../controls/select/select-old';
import { useSelector, useDispatch } from 'react-redux';
import { actions as associationActions } from '../../../services/store/association';
import { actions as modalActions } from '../../../services/store/modal';
import { getJSON } from '../../../services/utils';

const actions = { ...associationActions, ...modalActions };
const { Group } = Form;

export default function UserForm() {
  const dispatch = useDispatch();
  const mergeState = async (state) =>
    await dispatch(actions.mergeAssociation({ main: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));
  const resetAssociation = (_) => dispatch(actions.resetAssociation());

  const {
    submitted,
    source,
    loadingData,
    exposureSignature,
    studyOptions,
    strategyOptions,
    cancerOptions,
    rsSetOptions,
    study,
    strategy,
    cancer,
    rsSet,
  } = useSelector((state) => state.association.main);

  // populate controls on inital render
  useEffect(() => {
    if (!studyOptions.length) populateControls();
  }, []);

  // populate form
  async function populateControls() {
    mergeState({ loadingData: true });

    try {
      const exposureSignature = await getJSON(
        'Others/json/Exploring-Exposure.json'
      );

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
      const rsSet = 'COSMIC_v3_Signatures_GRCh37_SBS96'; // default

      const cancerOptions = [
        ...new Set(
          exposureSignature
            .filter((data) => data.Study == study && data.Dataset == strategy)
            .map((data) => data.Cancer_Type)
        ),
      ].sort();
      const cancer = 'Lung-AdenoCA'; // default

      mergeState({
        exposureSignature,
        study,
        studyOptions,
        strategy,
        strategyOptions,
        cancer,
        cancerOptions,
        rsSet,
        rsSetOptions,
      });
    } catch (error) {
      mergeError(error);
    }

    mergeState({ loadingData: false });
  }

  function handleReset() {
    const params = {
      source,
      exposureSignature,
      studyOptions,
      strategyOptions,
      cancerOptions,
      rsSetOptions,
      study: 'PCAWG',
      strategy: 'WGS',
      rsSet: 'COSMIC_v3_Signatures_GRCh37_SBS96',
      cancer: 'Lung-AdenoCA',
    };
    resetAssociation();
    mergeState(params);
  }

  async function handleLoadData() {
    const loadData = async () =>
      (
        await fetch(`web/associationWrapper`, {
          method: 'POST',
          headers: {
            Accept: 'application/json',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            fn: 'loadData',
            args: { study, strategy, rsSet, cancer },
          }),
        })
      ).json();

    mergeState({ loadingData: true });
    try {
      const [assocVarData, exposureVariableData] = await Promise.all([
        getJSON(`Association/PCAWG_vardata.json`),
        loadData(),
      ]);

      const { expVarList, error } = exposureVariableData;
      if (error) throw error.message;

      mergeState({
        assocVarData,
        expVarList,
        submitted: true,
      });
    } catch (error) {
      mergeError(error);
    }
    mergeState({ loadingData: false });
  }

  function handleStudy(study) {
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
    const rsSet = rsSetOptions[0];

    const cancerOptions = [
      ...new Set(
        exposureSignature
          .filter((data) => data.Study == study && data.Dataset == strategy)
          .map((data) => data.Cancer_Type)
      ),
    ];

    handleSet(rsSet);

    mergeState({
      study: study,
      strategy: strategy,
      strategyOptions: strategyOptions,
      cancer: cancerOptions[0],
      cancerOptions: cancerOptions,
      rsSetOptions: rsSetOptions,
      rsSet: rsSet,
    });
  }

  function handleStrategy(strategy) {
    const rsSetOptions = [
      ...new Set(
        exposureSignature
          .filter((row) => row.Study == study && row.Dataset == strategy)
          .map((row) => row.Signature_set_name)
      ),
    ];
    const rsSet = rsSetOptions[0];

    const cancerOptions = [
      ...new Set(
        exposureSignature
          .filter((data) => data.Study == study && data.Dataset == strategy)

          .map((data) => data.Cancer_Type)
      ),
    ].sort();

    handleSet(rsSet);

    mergeState({
      strategy: strategy,
      cancer: cancerOptions[0],
      cancerOptions: cancerOptions,
      rsSetOptions: rsSetOptions,
      rsSet: rsSet,
    });
  }

  function handleSet(set) {
    const rsSetOptions = [
      ...new Set(
        exposureSignature
          .filter((row) => row.Signature_set_name == set)
          .map((row) => row.Signature_name)
      ),
    ];

    mergeState({
      rsSet: set,
      rsSetOptions: rsSetOptions,
    });
  }

  return (
    <Form>
      <Row>
        <Col>
          <Group>
            <CustomSelect
              disabled={loadingData || submitted}
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
            <CustomSelect
              disabled={loadingData || submitted}
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
            <CustomSelect
              disabled={loadingData || submitted}
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
            <CustomSelect
              className="mb-4"
              disabled={loadingData || submitted}
              id="prevalenceCancerType"
              label="Cancer Type or Group"
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
        <Col md="6">
          <Button
            disabled={loadingData}
            className="w-100 mb-3"
            variant="secondary"
            onClick={() => handleReset()}
          >
            Reset
          </Button>
        </Col>
        <Col md="6">
          <Button
            disabled={loadingData || submitted}
            className="w-100"
            variant="primary"
            onClick={() => handleLoadData()}
          >
            Load Data
          </Button>
        </Col>
      </Row>
    </Form>
  );
}
