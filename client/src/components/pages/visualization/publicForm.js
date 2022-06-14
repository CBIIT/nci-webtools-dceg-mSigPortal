import React, { useEffect } from 'react';
import { Form, Button, Row, Col, Popover } from 'react-bootstrap';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import CustomSelect from '../../controls/select/select-old';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../services/store/visualization';
import { actions as modalActions } from '../../../services/store/modal';
import axios from 'axios';

const actions = { ...visualizationActions, ...modalActions };

export default function PublicForm() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.visualization);

  const mergeState = (state) =>
    dispatch(actions.mergeVisualization({ state: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));
  const resetVisualization = (_) => dispatch(actions.resetVisualization());

  const { study, cancer, strategy, data, loading } = store.publicForm;
  const { source, submitted } = store.state;

  useEffect(() => {
    if (!data.length && !loading) getPublicDataOptions();
  }, [data, source]);

  async function handleSubmit() {
    const args = {
      study,
      cancer,
      strategy,
    };

    mergeState({
      loading: {
        active: true,
        // content: 'Loading Public Data',
        // showIndicator: true,
      },
    });
    try {
      const response = await fetch(`web/visualizationWrapper`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ fn: 'getPublicData', args }),
      });

      if (response.ok) {
        const { output, projectID } = await response.json();
        if (output.error) throw output.error;
        if (output.uncaughtError) throw output.uncaughtError;

        mergeState({
          svgList: output,
          projectID: projectID,
        });
      } else if (response.status == 504) {
        mergeState({
          error: 'Please Reset Your Parameters and Try again.',
        });
        mergeError({
          visible: true,
          message:
            'Your submission has timed out. Please try again by submitting this job to a queue instead.',
        });
      } else {
        mergeState({
          error: 'Please Reset Your Parameters and Try again.',
        });
      }
    } catch (err) {
      mergeError(err.message);
      mergeState({
        error: 'Please Reset Your Parameters and Try again.',
      });
    }
    mergeState({ loading: { active: false } });
  }

  function handleReset() {
    const params = {
      data: data,
      // studyOptions: studyOptions,
      // study: studyOptions[0],
      // cancerOptions: cancerOptions,
      // cancer: 'Lung-AdenoCA',
      // pubExperimentOptions: pubExperimentOptions,
      // strategy: pubExperimentOptions[0],
    };
    window.location.hash = '#/visualization';
    resetVisualization();
    mergeState(params);
  }

  const studyOptions = [...new Set(data.map((data) => data.Study))];
  // default study

  const cancerOptions = [
    ...new Set(
      data.filter((data) => data.Study == study).map((data) => data.Cancer_Type)
    ),
  ];
  //  default cancer type

  const pubExperimentOptions = [
    ...new Set(
      data
        .filter((data) => data.Study == study && data.Cancer_Type == cancer)
        .map((data) => data.Dataset)
    ),
  ];

  async function getPublicDataOptions() {
    mergeState({ loading: true });
    try {
      const { data } = await axios.post('web/getFileS3', {
        path: `Others/json/Visualization-Public.json`,
      });

      mergeState({ data: data });

      // mergeCosineSimilarity({
      //   pubStudy: study,
      //   pubcancer: cancer,
      //   pubcancerOptions: cancerOptions,
      // });

      // mergeProfileComparison({
      //   pubStudy: study,
      //   pubcancer: cancer,
      //   pubcancerOptions: cancerOptions,
      // });

      // mergePCA({
      //   pubStudy: study,
      //   pubcancer: cancer,
      //   pubcancerOptions: cancerOptions,
      // });
    } catch (err) {
      mergeError(err.message);
    }
    mergeState({ loading: false });
  }

  function handleStudyChange(study) {
    const cancerOptions = [
      ...new Set(
        data
          .filter((data) => data.Study == study)
          .map((data) => data.Cancer_Type)
      ),
    ];

    const esOptions = [
      ...new Set(
        data
          .filter(
            (data) =>
              data.Study == study && data.Cancer_Type == cancerOptions[0]
          )
          .map((data) => data.Dataset)
      ),
    ];

    mergeState({
      study: study,
      cancer: cancerOptions[0],
      cancerOptions: cancerOptions,
      strategy: esOptions[0],
      pubExperimentOptions: esOptions,
    });
  }

  function handleCancerChange(cancer) {
    const pubExperimentOptions = [
      ...new Set(
        data
          .filter((data) => data.Study == study && data.Cancer_Type == cancer)
          .map((data) => data.Dataset)
      ),
    ];

    mergeState({
      cancer: cancer,
      strategy: pubExperimentOptions[0],
      pubExperimentOptions: pubExperimentOptions,
    });
  }

  return (
    <Form>
      <LoadingOverlay active={false} />
      <CustomSelect
        className="mb-2"
        id="publicFormStudy"
        label="Study"
        disabled={submitted}
        value={study}
        options={studyOptions}
        onChange={handleStudyChange}
      />
      <CustomSelect
        className="mb-2"
        id="publicFromcancer"
        label="Cancer Type or Group"
        disabled={submitted}
        value={cancer}
        options={cancerOptions}
        onChange={handleCancerChange}
      />
      <CustomSelect
        className="mb-4"
        id="publicFormStrategy"
        label="Experimental Strategy"
        disabled={submitted}
        value={strategy}
        options={pubExperimentOptions}
        onChange={(strategy) =>
          mergeState({
            strategy: strategy,
          })
        }
      />

      <Row>
        <Col>
          <Button
            disabled={loading.active || loading}
            className="w-100"
            variant="secondary"
            onClick={() => handleReset()}
          >
            Reset
          </Button>
        </Col>
        <Col>
          <Button
            disabled={loading.active || loading || submitted}
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
