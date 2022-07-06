import { Form, Button, Row, Col } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Select from '../../../controls/select/selectForm';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../../services/store/visualization';
import { actions as modalActions } from '../../../../services/store/modal';
import {
  useSubmitWebPublicMutation,
  useGetPublicDataOptionsQuery,
} from './apiSlice';

const actions = { ...visualizationActions, ...modalActions };

export default function PublicForm() {
  const store = useSelector((state) => state.visualization);
  const { submitted } = store.main;

  const dispatch = useDispatch();
  const mergeState = (state) =>
    dispatch(actions.mergeVisualization({ publicForm: state }));
  const resetVisualization = (_) => dispatch(actions.resetVisualization());
  const mergeMain = (state) =>
    dispatch(actions.mergeVisualization({ main: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const { data, error, isFetching } = useGetPublicDataOptionsQuery();
  const [handleSubmitWeb, { isLoading, reset: resetWeb }] =
    useSubmitWebPublicMutation();

  const defaultValues = {
    study: { label: 'PCAWG', value: 'PCAWG' },
    cancer: { label: 'Lung-AdenoCA', value: 'Lung-AdenoCA' },
    strategy: { label: 'WGS', value: 'WGS' },
  };

  const {
    control,
    handleSubmit,
    reset: resetForm,
    setValue,
    watch,
  } = useForm({ defaultValues });

  const formStudy = watch('study');
  const formCancer = watch('cancer');

  function handleReset() {
    window.location.hash = '#/visualization';
    resetForm();
    resetWeb();
    resetVisualization();
  }

  async function onSubmit(data) {
    try {
      mergeMain({ submitted: true, loading: { active: true } });
      mergeState(data);
      const params = {
        study: data.study.value,
        cancer: data.cancer.value,
        strategy: data.strategy.value,
      };
      const { data: samples, projectID } = await handleSubmitWeb(params).unwrap();

      mergeMain({ samples, projectID: projectID });
    } catch (error) {
      if (error.originalStatus == 504) {
        mergeMain({
          error: 'Please Reset Your Parameters and Try again.',
        });
        mergeError({
          visible: true,
          message:
            'Your submission has timed out. Please try again by submitting this job to a queue instead.',
        });
      } else {
        mergeMain({
          error: 'Please Reset Your Parameters and Try again.',
        });
        mergeError(error.data);
      }
    }
    mergeMain({ loading: { active: false } });
  }

  const studyOptions = data
    ? Object.keys(data).map((e) => ({ label: e, value: e }))
    : [];

  const cancerOptions = (study) => {
    if (data && study.value) {
      const options = Object.keys(data[study.value]).map((e) => ({
        label: e,
        value: e,
      }));
      return options;
    } else {
      return [];
    }
  };

  const strategyOptions = (study, cancer) => {
    if (data && study.value && cancer.value) {
      const strategies = data[study.value][cancer.value];
      if (strategies) {
        const options = Object.values(strategies).map((e) => ({
          label: e,
          value: e,
        }));
        return options;
      } else {
        return [];
      }
    } else {
      return [];
    }
  };

  function handleStudyChange(study) {
    const cancers = cancerOptions(study);
    const strategies = strategyOptions(study, cancers[0]);

    setValue('study', study);
    setValue('cancer', cancers[0]);
    setValue('strategy', strategies[0]);
  }

  function handleCancerChange(cancer) {
    const strategies = strategyOptions(formStudy, cancer);

    setValue('cancer', cancer);
    setValue('strategy', strategies[0]);
  }

  return (
    <Form onSubmit={handleSubmit(onSubmit)}>
      {error && <p>There was an error retrieving public data options</p>}
      <Select
        className="mb-2"
        name="study"
        label="Study"
        disabled={submitted}
        options={studyOptions}
        control={control}
        onChange={handleStudyChange}
      />
      <Select
        className="mb-2"
        name="cancer"
        label="Cancer Type or Group"
        disabled={submitted}
        options={cancerOptions(formStudy)}
        control={control}
        onChange={handleCancerChange}
      />
      <Select
        className="mb-4"
        name="strategy"
        label="Experimental Strategy"
        disabled={submitted}
        options={strategyOptions(formStudy, formCancer)}
        control={control}
      />

      <Row>
        <Col>
          <Button
            disabled={isFetching || isLoading}
            className="w-100"
            variant="secondary"
            onClick={() => handleReset()}
          >
            Reset
          </Button>
        </Col>
        <Col>
          <Button
            disabled={isFetching || isLoading || submitted}
            className="w-100"
            variant="primary"
            type="submit"
          >
            Submit
          </Button>
        </Col>
      </Row>
    </Form>
  );
}
