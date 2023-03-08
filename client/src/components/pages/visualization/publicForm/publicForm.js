import { Form, Button, Row, Col } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import Select from '../../../controls/select/selectForm';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../../services/store/visualization';
import { actions as modalActions } from '../../../../services/store/modal';
import {
  resetVisualizationApi,
  useSeqmatrixOptionsQuery,
} from '../../../../services/store/rootApi';
import { usePublicMatrixMutation } from './apiSlice';

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

  const { data, error, isFetching } = useSeqmatrixOptionsQuery();
  const [fetchMatrix, { isLoading }] = usePublicMatrixMutation();

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
    dispatch(resetVisualizationApi);
    resetVisualization();
  }

  async function* paginateQuery(endpoint, params) {
    const limit = 1000000;
    let offset = 0;
    let result = [];

    do {
      result = await endpoint({ ...params, limit, offset }).unwrap();
      offset += limit;
      yield result;
    } while (result.length >= limit);
    return result;
  }

  async function onSubmit(data) {
    try {
      const cancers = cancerOptions(data.study);
      mergeMain({ submitted: true, loading: { active: true } });
      mergeState({...data, cancers });
      const params = {
        study: data.study.value,
        cancer: data.cancer.value,
        strategy: data.strategy.value,
      };

      // let matrixData = [];
      // for await (const data of paginateQuery(fetchMatrix, params)) {
      //   matrixData = [...matrixData, ...data];
      // }
      const matrixData = await fetchMatrix(params).unwrap();

      mergeMain({ matrixData, id: crypto.randomUUID() });
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
    ? [...new Set(data.map((e) => e.study))].sort().map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  const cancerOptions = (study) => {
    if (data && study.value) {
      const options = [
        ...[
          ...new Set(
            data.filter((e) => e.study == study.value).map((e) => e.cancer)
          ),
        ]
          .sort()
          .map((e) => ({
            label: e,
            value: e,
          })),
      ];
      if (study.value == 'PCAWG' || study.value == 'TCGA') {
        return [
          ...options,
          {
            label: 'PanCancer',
            value: '*ALL',
          },
        ];
      } else {
        return options;
      }
    } else {
      return [];
    }
  };

  const strategyOptions = (study, cancer) => {
    if (data && study.value) {
      const options =
        cancer.value == '*ALL'
          ? [
              ...new Set(
                data
                  .filter((e) => e.study == study.value)
                  .map((e) => e.strategy)
              ),
            ]
          : [
              ...new Set(
                data
                  .filter(
                    (e) => e.study == study.value && e.cancer == cancer.value
                  )
                  .map((e) => e.strategy)
              ),
            ];
      return options.sort().map((e) => ({
        label: e,
        value: e,
      }));
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
      <LoadingOverlay active={isFetching} />
      {error && <p>There was an error retrieving public data options</p>}
      <Select
        className="mb-2"
        name="study"
        label="Study"
        disabled={submitted || isFetching}
        options={studyOptions}
        control={control}
        onChange={handleStudyChange}
      />
      <Select
        className="mb-2"
        name="cancer"
        label="Cancer Type or Group"
        disabled={submitted || isFetching}
        options={cancerOptions(formStudy)}
        control={control}
        onChange={handleCancerChange}
      />
      <Select
        className="mb-4"
        name="strategy"
        label="Experimental Strategy"
        disabled={submitted || isFetching}
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
