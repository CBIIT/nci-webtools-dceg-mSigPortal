import React, { useEffect } from 'react';
import { Form, Button, Row, Col, Popover } from 'react-bootstrap';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Select from '../../controls/select/select';
import { useSelector, useDispatch } from 'react-redux';
import {
  actions as visualizationActions,
  getInitialState,
} from '../../../services/store/visualization';
import { actions as modalActions } from '../../../services/store/modal';

const actions = { ...visualizationActions, ...modalActions };

export default function PublicForm() {
  const dispatch = useDispatch();
  const visualization = useSelector((state) => state.visualization);

  const mergeVisualize = (state) =>
    dispatch(actions.mergeVisualization({ visualize: state }));
  const mergeResults = (state) =>
    dispatch(actions.mergeVisualization({ results: state }));
  const mergeProfilerSummary = (state) =>
    dispatch(actions.mergeVisualization({ profilerSummary: state }));
  const mergeMutationalProfiles = (state) =>
    dispatch(actions.mergeVisualization({ mutationalProfiles: state }));
  const mergeMutationalPattern = (state) =>
    dispatch(actions.mergeVisualization({ mutationalPattern: state }));
  const mergeCosineSimilarity = (state) =>
    dispatch(actions.mergeVisualization({ cosineSimilarity: state }));
  const mergeProfileComparison = (state) =>
    dispatch(actions.mergeVisualization({ profileComparison: state }));
  const mergePCA = (state) =>
    dispatch(actions.mergeVisualization({ pca: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));
  const resetVisualization = (_) => dispatch(actions.resetVisualization());

  const {
    study,
    studyOptions,
    cancerType,
    cancerTypeOptions,
    pubExperimentalStrategy,
    pubExperimentOptions,
    pDataOptions,
    submitted,
    loading,
    loadingPublic,
    source,
  } = visualization.visualize;

  useEffect(() => {
    if (!pDataOptions.length && !loadingPublic) getPublicDataOptions();
  }, [pDataOptions, source]);

  async function handleSubmit() {
    // disable parameters after submit
    mergeVisualize({ visualize: { submitted: true } });

    const args = {
      study: study,
      cancerType: cancerType,
      experimentalStrategy: pubExperimentalStrategy,
    };

    mergeVisualize({
      loading: {
        active: true,
        content: 'Loading Public Data',
        showIndicator: true,
      },
    });
    try {
      const response = await fetch(`api/getPublicData`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(args),
      });

      if (response.ok) {
        const { svgList, projectID } = await response.json();

        mergeResults({
          svgList: svgList,
          projectID: projectID,
        });
      } else if (response.status == 504) {
        mergeResults({
          error: 'Please Reset Your Parameters and Try again.',
        });
        mergeError({
          visible: true,
          message:
            'Your submission has timed out. Please try again by submitting this job to a queue instead.',
        });
      } else {
        mergeResults({
          error: 'Please Reset Your Parameters and Try again.',
        });
      }
    } catch (err) {
      mergeError(err.message);
      mergeResults({
        error: 'Please Reset Your Parameters and Try again.',
      });
    }
    mergeVisualize({ loading: { active: false } });
  }

  function handleReset() {
    const params = {
      pDataOptions: pDataOptions,
      studyOptions: studyOptions,
      study: studyOptions[0],
      cancerTypeOptions: cancerTypeOptions,
      cancerType: 'Lung-AdenoCA',
      pubExperimentOptions: pubExperimentOptions,
      pubExperimentalStrategy: pubExperimentOptions[0],
    };
    resetVisualization();
    mergeVisualize(params);
  }

  async function getPublicDataOptions() {
    mergeVisualize({ loadingPublic: true });
    try {
      const pDataOptions = await (
        await fetch(`api/public/Others/json/Visualization-Public.json`)
      ).json();

      const studyOptions = [...new Set(pDataOptions.map((data) => data.Study))];
      // default study
      const study = 'PCAWG';

      const cancerTypeOptions = [
        ...new Set(
          pDataOptions
            .filter((data) => data.Study == study)
            .map((data) => data.Cancer_Type)
        ),
      ];
      //  default cancer type
      const cancer = 'Lung-AdenoCA';

      const pubExperimentOptions = [
        ...new Set(
          pDataOptions
            .filter((data) => data.Study == study && data.Cancer_Type == cancer)
            .map((data) => data.Dataset)
        ),
      ];

      mergeVisualize({
        pDataOptions: pDataOptions,
        study: study,
        studyOptions: studyOptions,
        cancerType: cancer,
        cancerTypeOptions: cancerTypeOptions,
        pubExperimentalStrategy: pubExperimentOptions[0],
        pubExperimentOptions: pubExperimentOptions,
      });

      mergeCosineSimilarity({
        pubStudy: study,
        pubCancerType: cancer,
        pubCancerTypeOptions: cancerTypeOptions,
      });

      mergeProfileComparison({
        pubStudy: study,
        pubCancerType: cancer,
        pubCancerTypeOptions: cancerTypeOptions,
      });

      mergePCA({
        pubStudy: study,
        pubCancerType: cancer,
        pubCancerTypeOptions: cancerTypeOptions,
      });
    } catch (err) {
      mergeError(err.message);
    }
    mergeVisualize({ loadingPublic: false });
  }

  function handleStudyChange(study) {
    const cancerTypeOptions = [
      ...new Set(
        pDataOptions
          .filter((data) => data.Study == study)
          .map((data) => data.Cancer_Type)
      ),
    ];

    const esOptions = [
      ...new Set(
        pDataOptions
          .filter(
            (data) =>
              data.Study == study && data.Cancer_Type == cancerTypeOptions[0]
          )
          .map((data) => data.Dataset)
      ),
    ];

    mergeVisualize({
      study: study,
      cancerType: cancerTypeOptions[0],
      cancerTypeOptions: cancerTypeOptions,
      pubExperimentalStrategy: esOptions[0],
      pubExperimentOptions: esOptions,
    });
  }

  function handleCancerChange(cancer) {
    const pubExperimentOptions = [
      ...new Set(
        pDataOptions
          .filter((data) => data.Study == study && data.Cancer_Type == cancer)
          .map((data) => data.Dataset)
      ),
    ];

    mergeVisualize({
      cancerType: cancer,
      pubExperimentalStrategy: pubExperimentOptions[0],
      pubExperimentOptions: pubExperimentOptions,
    });
  }

  return (
    <Form>
      <LoadingOverlay active={loadingPublic} />
      <Select
        className="mb-2"
        id="publicFormStudy"
        label="Study"
        disabled={submitted}
        value={study}
        options={studyOptions}
        onChange={handleStudyChange}
      />
      <Select
        className="mb-2"
        id="publicFromCancerType"
        label="Cancer Type"
        disabled={submitted}
        value={cancerType}
        options={cancerTypeOptions}
        onChange={handleCancerChange}
      />
      <Select
        className="mb-4"
        id="publicFormStrategy"
        label="Experimental Strategy"
        disabled={submitted}
        value={pubExperimentalStrategy}
        options={pubExperimentOptions}
        onChange={(pubExperimentalStrategy) =>
          mergeVisualize({
            pubExperimentalStrategy: pubExperimentalStrategy,
          })
        }
      />

      <Row>
        <Col lg="6">
          <Button
            disabled={loading.active || loadingPublic}
            className="w-100"
            variant="secondary"
            onClick={(e) => handleReset()}
          >
            Reset
          </Button>
        </Col>
        <Col lg="6">
          <Button
            disabled={loading.active || loadingPublic || submitted}
            className="w-100"
            variant="primary"
            type="button"
            onClick={(e) => handleSubmit(e)}
          >
            Submit
          </Button>
        </Col>
      </Row>
    </Form>
  );
}
