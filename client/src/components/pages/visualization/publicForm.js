import React, { useEffect } from 'react';
import { Form, Button, Row, Col } from 'react-bootstrap';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Select from '../../controls/select/selectForm';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../services/store/visualization';
import { actions as modalActions } from '../../../services/store/modal';
import axios from 'axios';

const actions = { ...visualizationActions, ...modalActions };

export default function PublicForm() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.visualization);

  const mergeState = (state) =>
    dispatch(actions.mergeVisualization({ publicForm: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));
  const resetVisualization = (_) => dispatch(actions.resetVisualization());

  const { study, cancer, strategy, data, loading, defaultOptions } =
    store.publicForm;
  const { source, submitted } = store.main;

  // get options data
  useEffect(() => {
    if (!data.length && !loading && source == 'public') getPublicDataOptions();
  }, [data, source]);

  // set default options
  useEffect(() => {
    if (data.length && study == null) mergeState({ ...defaultOptions });
  }, [data]);

  const studyOptions = () =>
    data.length
      ? [...new Set(data.map((e) => e.Study))].map((e) => ({
          label: e,
          value: e,
        }))
      : [];

  const cancerOptions = (study) =>
    study && data.length
      ? [
          ...new Set(
            data.filter((e) => e.Study == study.value).map((e) => e.Cancer_Type)
          ),
        ].map((e) => ({ label: e, value: e }))
      : [];

  const strategyOptions = (study, cancer) =>
    study && cancer && data.length
      ? [
          ...new Set(
            data
              .filter(
                (e) => e.Study == study.value && e.Cancer_Type == cancer.value
              )
              .map((e) => e.Dataset)
          ),
        ].map((e) => ({ label: e, value: e }))
      : [];

  async function handleSubmit() {
    dispatch(actions.mergeVisualization({ main: { submitted: true } }));
  }

  function handleReset() {
    const cacheData = { data: [...data], ...defaultOptions };
    window.location.hash = '#/visualization';
    resetVisualization();
    mergeState(cacheData);
  }

  async function getPublicDataOptions() {
    mergeState({ loading: true });
    try {
      const { data } = await axios.post('web/getFileS3', {
        path: `Others/json/Visualization-Public.json`,
      });

      mergeState({ data: data });
    } catch (err) {
      mergeError(err.message);
    }
    mergeState({ loading: false });
  }

  function handleStudyChange(study) {
    const cancers = cancerOptions(study);
    const strategies = strategyOptions(study, cancers[0]);

    mergeState({
      study: study,
      cancer: cancers[0],
      strategy: strategies[0],
    });
  }

  function handleCancerChange(cancer) {
    const strategies = strategyOptions(study, cancer);

    mergeState({
      cancer: cancer,
      strategy: strategies[0],
    });
  }

  function handleStrategyChange(strategy) {
    mergeState({ strategy });
  }

  return (
    <Form>
      <LoadingOverlay active={loading} />
      <Select
        className="mb-2"
        id="publicStudy"
        label="Study"
        disabled={submitted}
        value={study}
        options={studyOptions()}
        onChange={handleStudyChange}
      />
      <Select
        className="mb-2"
        id="publicCancer"
        label="Cancer Type or Group"
        disabled={submitted}
        value={cancer}
        options={cancerOptions(study)}
        onChange={handleCancerChange}
      />
      <Select
        className="mb-4"
        id="publicStrategy"
        label="Experimental Strategy"
        disabled={submitted}
        value={strategy}
        options={strategyOptions(study, cancer)}
        onChange={handleStrategyChange}
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
            onClick={() => handleSubmit()}
          >
            Submit
          </Button>
        </Col>
      </Row>
    </Form>
  );
}
