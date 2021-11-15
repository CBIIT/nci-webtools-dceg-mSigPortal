import React, { useEffect, useState } from 'react';
import { Form, Row, Col, Button, Modal } from 'react-bootstrap';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Select from '../../controls/select/select';
import { useSelector, useDispatch } from 'react-redux';
import { actions as associationActions } from '../../../services/store/association';
import { actions as modalActions } from '../../../services/store/modal';
import { getJSON } from '../../../services/utils';

const actions = { ...associationActions, ...modalActions };
const { Group } = Form;

export default function PublicForm() {
  const dispatch = useDispatch();
  const mergeState = async (state) =>
    dispatch(actions.mergeAssociation({ associationState: state }));
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
  } = useSelector((state) => state.association.associationState);

  const [missingModal, setMissing] = useState(false);

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
      ];
      const cancer = 'Panc-AdenoCA'; // default

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
    // tempoary - display warning for missing or incomplete studies
    if (study == 'Mutographs-ESCC') {
      setMissing(true);
      mergeState({ submitted: true });
    } else {
      const getAssocVarData = async () =>
        (
          await fetch(`api/associationWrapper`, {
            method: 'POST',
            headers: {
              Accept: 'application/json',
              'Content-Type': 'application/json',
            },
            body: JSON.stringify({
              fn: 'getAssocVarData',
              args: { study, cancer },
            }),
          })
        ).json();
      const getExpVarData = async () =>
        (
          await fetch(`api/associationWrapper`, {
            method: 'POST',
            headers: {
              Accept: 'application/json',
              'Content-Type': 'application/json',
            },
            body: JSON.stringify({
              fn: 'getExpVarData',
              args: { study, strategy, rsSet, cancer },
            }),
          })
        ).json();

      mergeState({ loadingData: true });
      try {
        const [assocResponse, expResponse] = await Promise.all([
          getAssocVarData(),
          getExpVarData(),
        ]);

        const {
          projectID,
          output: assocOutput,
          uncaughtError: assocError,
        } = assocResponse;
        const { output: expOutput, uncaughtError: expError } = expResponse;

        if (assocError) throw assocError;
        if (expError) throw expError;

        mergeState({
          submitted: true,
          assocVarData: assocOutput.assocVarData,
          assocFullDataPath: assocOutput.fullDataPath,
          expVarList: expOutput.expVarList,
          displayTab: 'univariable',
          openSidebar: false,
        });
        dispatch(actions.mergeAssociation({ univariable: { projectID } }));
      } catch (error) {
        mergeError(error);
      }
      mergeState({ loadingData: false });
    }
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
    ];

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
    mergeState({ rsSet: set });
  }

  return (
    <>
      <Modal
        data-testid="missingModal"
        show={missingModal}
        onHide={() => setMissing(false)}
      >
        <Modal.Body style={{ backgroundColor: '#fafafa' }}>
          <p>Selected data is currently unavailable.</p>
          <p>
            {`${study}-${rsSet}-${cancer} will be added in a future release.`}
          </p>
        </Modal.Body>

        <Modal.Footer
          className="d-flex justify-content-center border-0"
          style={{ backgroundColor: '#fafafa' }}
        >
          <Button variant="outline-secondary" onClick={() => setMissing(false)}>
            Close
          </Button>
        </Modal.Footer>
      </Modal>

      <Form>
        <LoadingOverlay active={!studyOptions.length} />
        <Row>
          <Col>
            <Group>
              <Select
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
              <Select
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
              <Select
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
              <Select
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
              Load Study
            </Button>
          </Col>
        </Row>
      </Form>
    </>
  );
}
