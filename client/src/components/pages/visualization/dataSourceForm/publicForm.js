import React from 'react';
import { Form, Button, Row, Col } from 'react-bootstrap';
import { useRecoilState, useRecoilValue, useResetRecoilState } from 'recoil';
// import { useForm, Controller } from 'react-hook-form';
import Select from 'react-select';
import {
  visualizationState,
  publicFormState,
  publicFormOptions,
} from '../visualization.state';

export default function PublicForm() {
  const [formState, setForm] = useRecoilState(publicFormState);
  const mergeForm = (state) => setForm({ ...formState, ...state });
  const reset = useResetRecoilState(publicFormState);

  const { study, cancer, strategy } = formState;
  const { studyOptions, cancerOptions, strategyOptions } =
    useRecoilValue(publicFormOptions);
  const { submitted } = useRecoilValue(visualizationState);

  function handleSelect(value, action) {
    console.log(action);
    const { name } = action;
    mergeForm({ [name]: value });
  }

  async function handleSubmit() {
    const args = {
      study: study,
      cancer: cancer,
      strategy,
    };

    // mergeState({
    //   loading: {
    //     active: true,
    //     // content: 'Loading Public Data',
    //     // showIndicator: true,
    //   },
    // });
    // try {
    //   const response = await fetch(`web/visualizationWrapper`, {
    //     method: 'POST',
    //     headers: {
    //       Accept: 'application/json',
    //       'Content-Type': 'application/json',
    //     },
    //     body: JSON.stringify({ fn: 'getPublicData', args }),
    //   });

    //   if (response.ok) {
    //     const { output, projectID } = await response.json();
    //     if (output.error) throw output.error;
    //     if (output.uncaughtError) throw output.uncaughtError;

    //     mergeState({
    //       svgList: output,
    //       projectID: projectID,
    //     });
    //   } else if (response.status == 504) {
    //     mergeState({
    //       error: 'Please Reset Your Parameters and Try again.',
    //     });
    //     mergeError({
    //       visible: true,
    //       message:
    //         'Your submission has timed out. Please try again by submitting this job to a queue instead.',
    //     });
    //   } else {
    //     mergeState({
    //       error: 'Please Reset Your Parameters and Try again.',
    //     });
    //   }
    // } catch (err) {
    //   mergeError(err.message);
    //   mergeState({
    //     error: 'Please Reset Your Parameters and Try again.',
    //   });
    // }
    // mergeState({ loading: { active: false } });
  }

  function handleReset() {
    // const params = {
    //   data: data,
    //   studyOptions: studyOptions,
    //   study: studyOptions[0],
    //   cancerOptions: cancerOptions,
    //   cancer: 'Lung-AdenoCA',
    //   pubExperimentOptions: pubExperimentOptions,
    //   strategy: pubExperimentOptions[0],
    // };
    reset();
    window.location.hash = '#/visualization';
    // resetVisualization();
    // mergeState(params);
  }

  // async function getPublicDataOptions() {
  //   try {
  //     // const data = await getJSON(
  //     //   `Others/json/Visualization-Public.json`
  //     // );
  //     const studyOptions = [...new Set(data.map((data) => data.Study))];
  //     // default study
  //     const study = 'PCAWG';
  //     const cancerOptions = [
  //       ...new Set(
  //         data
  //           .filter((data) => data.Study == study)
  //           .map((data) => data.Cancer_Type)
  //       ),
  //     ];
  //     //  default cancer type
  //     const cancer = 'Lung-AdenoCA';
  //     const pubExperimentOptions = [
  //       ...new Set(
  //         data
  //           .filter((data) => data.Study == study && data.Cancer_Type == cancer)
  //           .map((data) => data.Dataset)
  //       ),
  //     ];
  //     mergeState({
  //       data: data,
  //       study: study,
  //       studyOptions: studyOptions,
  //       cancer: cancer,
  //       cancerOptions: cancerOptions,
  //       strategy: pubExperimentOptions[0],
  //       pubExperimentOptions: pubExperimentOptions,
  //     });
  //     mergeCosineSimilarity({
  //       pubStudy: study,
  //       pubcancer: cancer,
  //       pubcancerOptions: cancerOptions,
  //     });
  //     mergeProfileComparison({
  //       pubStudy: study,
  //       pubcancer: cancer,
  //       pubcancerOptions: cancerOptions,
  //     });
  //     mergePCA({
  //       pubStudy: study,
  //       pubcancer: cancer,
  //       pubcancerOptions: cancerOptions,
  //     });
  //   } catch (err) {
  //     mergeError(err.message);
  //   }
  // }

  function selectValue(value, options) {
    if (options.some((e) => e.value == value.value)) {
      return value;
    } else {
      return options[0];
    }
  }

  return (
    <Form>
      <Form.Group controlId="study">
        <Form.Label>Study</Form.Label>
        <Select
          className="mb-2"
          inputId="study"
          name="study"
          disabled={submitted}
          value={selectValue(study, studyOptions)}
          options={studyOptions}
          onChange={handleSelect}
        />
      </Form.Group>
      <Form.Group controlId="cancer">
        <Form.Label>Cancer Type or Group</Form.Label>
        <Select
          className="mb-2"
          inputId="cancer"
          name="cancer"
          disabled={submitted}
          value={selectValue(cancer, cancerOptions)}
          options={cancerOptions}
          onChange={handleSelect}
        />
      </Form.Group>
      <Form.Group controlId="strategy">
        <Form.Label>Experimental Strategy</Form.Label>
        <Select
          className="mb-4"
          inputId="strategy"
          name="strategy"
          disabled={submitted}
          value={selectValue(strategy, strategyOptions)}
          options={strategyOptions}
          onChange={handleSelect}
        />
      </Form.Group>

      <Row>
        <Col>
          <Button
            // disabled={loading.active || loadingPublic}
            className="w-100"
            variant="secondary"
            onClick={() => handleReset()}
          >
            Reset
          </Button>
        </Col>
        <Col>
          <Button
            disabled={submitted}
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
