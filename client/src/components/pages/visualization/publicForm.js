import React, { useEffect } from 'react';
import { Form, Button, Row, Col, Popover } from 'react-bootstrap';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { useSelector } from 'react-redux';
import './visualization.scss';
import {
  dispatchVisualize,
  dispatchVisualizeResults,
  dispatchError,
  getInitialState,
  dispatchMutationalProfiles,
  dispatchCosineSimilarity,
  dispatchProfileComparison,
  dispatchPCA,
  dispatchMutationalPattern,
  dispatchProfilerSummary,
} from '../../../services/store';
import Select from '../../controls/select/select';

export default function PublicForm() {
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
  } = useSelector((state) => state.visualize);
  const rootURL = window.location.pathname;

  useEffect(() => {
    if (!pDataOptions.length && !loadingPublic) getPublicDataOptions();
  }, [pDataOptions, source]);

  async function handleSubmit() {
    // disable parameters after submit
    dispatchVisualize({ submitted: true });

    const args = {
      study: study,
      cancerType: cancerType,
      experimentalStrategy: pubExperimentalStrategy,
    };

    dispatchVisualize({
      loading: {
        active: true,
        content: 'Loading Public Data',
        showIndicator: true,
      },
    });
    try {
      const response = await fetch(`${rootURL}getPublicData`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(args),
      });

      if (response.ok) {
        const { svgList, projectID } = await response.json();

        dispatchVisualizeResults({
          svgList: svgList,
          projectID: projectID,
        });
      } else if (response.status == 504) {
        dispatchVisualizeResults({
          error: 'Please Reset Your Parameters and Try again.',
        });
        dispatchError(
          'Your submission has timed out. Please try again by submitting this job to a queue instead.'
        );
      } else {
        dispatchVisualizeResults({
          error: 'Please Reset Your Parameters and Try again.',
        });
      }
    } catch (err) {
      dispatchError(err);
      dispatchVisualizeResults({
        error: 'Please Reset Your Parameters and Try again.',
      });
    }
    dispatchVisualize({ loading: { active: false } });
  }

  function handleReset() {
    const initialState = getInitialState();
    dispatchVisualize({
      ...initialState.visualize,
      pDataOptions: pDataOptions,
      studyOptions: studyOptions,
      study: studyOptions[0],
      cancerTypeOptions: cancerTypeOptions,
      cancerType: cancerTypeOptions[0],
      pubExperimentOptions: pubExperimentOptions,
      pubExperimentalStrategy: pubExperimentOptions[0],
    });
    dispatchProfilerSummary(initialState.profilerSummary);
    dispatchVisualizeResults(initialState.visualizeResults);
    dispatchMutationalProfiles(initialState.mutationalProfiles);
    dispatchMutationalPattern(initialState.mutationalPattern);
    dispatchCosineSimilarity(initialState.cosineSimilarity);
    dispatchProfileComparison(initialState.profileComparison);
    dispatchPCA(initialState.pca);
  }

  async function getPublicDataOptions() {
    dispatchVisualize({ loadingPublic: true });
    try {
      const response = await fetch(`${rootURL}getPublicDataOptions`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
      });
      if (!response.ok) {
        // console.log(await response.json());
      } else {
        const pDataOptions = await response.json();
        const studyOptions = [
          ...new Set(pDataOptions.map((data) => data.Study)),
        ];
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
              .filter(
                (data) => data.Study == study && data.Cancer_Type == cancer
              )
              .map((data) => data.Dataset)
          ),
        ];

        dispatchVisualize({
          pDataOptions: pDataOptions,
          study: study,
          studyOptions: studyOptions,
          cancerType: cancer,
          cancerTypeOptions: cancerTypeOptions,
          pubExperimentalStrategy: pubExperimentOptions[0],
          pubExperimentOptions: pubExperimentOptions,
        });

        dispatchCosineSimilarity({
          pubStudy: study,
          pubCancerType: cancer,
          pubCancerTypeOptions: cancerTypeOptions,
        });

        dispatchProfileComparison({
          pubStudy: study,
          pubCancerType: cancer,
          pubCancerTypeOptions: cancerTypeOptions,
        });

        dispatchPCA({
          pubStudy: study,
          pubCancerType: cancer,
          pubCancerTypeOptions: cancerTypeOptions,
        });
      }
    } catch (err) {
      dispatchError(err);
    }
    dispatchVisualize({ loadingPublic: false });
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

    dispatchVisualize({
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

    dispatchVisualize({
      cancerType: cancer,
      pubExperimentalStrategy: pubExperimentOptions[0],
      pubExperimentOptions: pubExperimentOptions,
    });
  }

  return (
    <Form>
      <LoadingOverlay active={loadingPublic} />
      <Select
        id="publicFormStudy"
        label="Study"
        disabled={submitted}
        value={study}
        options={studyOptions}
        onChange={handleStudyChange}
      />
      <Select
        id="publicFromCancerType"
        label="Cancer Type"
        disabled={submitted}
        value={cancerType}
        options={cancerTypeOptions}
        onChange={handleCancerChange}
      />
      <Select
        id="publicFormStrategy"
        label="Experimental Strategy"
        disabled={submitted}
        value={pubExperimentalStrategy}
        options={pubExperimentOptions}
        onChange={(pubExperimentalStrategy) =>
          dispatchVisualize({
            pubExperimentalStrategy: pubExperimentalStrategy,
          })
        }
      />

      <Row>
        <Col sm="6">
          <Button
            disabled={loading.active || loadingPublic}
            className="w-100"
            variant="secondary"
            onClick={(e) => handleReset()}
          >
            Reset
          </Button>
        </Col>
        <Col sm="6">
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
