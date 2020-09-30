import React, { useEffect } from 'react';
import { Form, Button, Row, Col, Popover } from 'react-bootstrap';
import Select from 'react-select';
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
const { Group, Label } = Form;

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
        content: 'Searching...',
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
    dispatchVisualize(initialState.visualize);
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

        const cancerTypeOptions = [
          ...new Set(
            pDataOptions
              .filter((data) => data.Study == studyOptions[0])
              .map((data) => data.Cancer_Type)
          ),
        ];
        const pubExperimentOptions = [
          ...new Set(
            pDataOptions
              .filter(
                (data) =>
                  data.Study == studyOptions[0] &&
                  data.Cancer_Type == cancerTypeOptions[0]
              )
              .map((data) => data.Dataset)
          ),
        ];

        dispatchVisualize({
          pDataOptions: pDataOptions,
          study: studyOptions[0],
          studyOptions: studyOptions,
          cancerType: cancerTypeOptions[0],
          cancerTypeOptions: cancerTypeOptions,
          pubExperimentalStrategy: pubExperimentOptions[0],
          pubExperimentOptions: pubExperimentOptions,
        });

        dispatchCosineSimilarity({
          pubStudy: studyOptions[0],
          pubCancerType: cancerTypeOptions[0],
          pubCancerTypeOptions: cancerTypeOptions,
        });

        dispatchProfileComparison({
          pubStudy: studyOptions[0],
          pubCancerType: cancerTypeOptions[0],
          pubCancerTypeOptions: cancerTypeOptions,
        });

        dispatchPCA({
          pubStudy: studyOptions[0],
          pubCancerType: cancerTypeOptions[0],
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
      <Group controlId="study">
        <Label>Study</Label>
        <Select
          inputId="study"
          isDisabled={submitted}
          options={studyOptions}
          value={[study]}
          onChange={(study) => handleStudyChange(study)}
          getOptionLabel={(option) => option}
          getOptionValue={(option) => option}
        />
      </Group>
      <Group controlId="cancerType">
        <Label>Cancer Type</Label>
        <Select
          inputId="cancerType"
          isDisabled={submitted}
          options={cancerTypeOptions}
          value={[cancerType]}
          onChange={(cancerType) => handleCancerChange(cancerType)}
          getOptionLabel={(option) => option}
          getOptionValue={(option) => option}
        />
      </Group>

      <Group controlId="strategy">
        <Label>Experimental Strategy</Label>
        <Select
          inputId="strategy"
          isDisabled={submitted}
          options={pubExperimentOptions}
          value={[pubExperimentalStrategy]}
          onChange={(pubExperimentalStrategy) =>
            dispatchVisualize({
              pubExperimentalStrategy: pubExperimentalStrategy,
            })
          }
          getOptionLabel={(option) => option}
          getOptionValue={(option) => option}
        />
      </Group>

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
