import { Form, Button, Row, Col } from 'react-bootstrap';
import { Controller, useForm } from 'react-hook-form';
import { LoadingOverlay } from '@/components/controls/loading-overlay/loading-overlay';
import { useSelector, useDispatch } from 'react-redux';
import { actions as exposureActions } from '@/services/store/exploration';
import { actions as modalActions } from '@/services/store/modal';
import {
  resetExplorationApi,
  useExposureOptionsQuery,
} from '@/services/store/rootApi';
import PublicDataFields from './publicDataFields';

const actions = { ...exposureActions, ...modalActions };

export default function PublicForm() {
  const store = useSelector((state) => state.exploration);
  const { submitted, openSidebar } = store.main;

  const dispatch = useDispatch();
  const mergeForm = (state) =>
    dispatch(actions.mergeExploration({ publicForm: state }));
  const mergeCompareForm = (state) =>
    dispatch(actions.mergeExploration({ comparePublicForm: state }));
  const resetExploration = (_) => dispatch(actions.resetExploration());
  const mergeMain = (state) =>
    dispatch(actions.mergeExploration({ main: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const {
    data: options,
    error: optionsError,
    isFetching,
  } = useExposureOptionsQuery();

  const publicDataDefaults = {
    study: { label: 'PCAWG', value: 'PCAWG' },
    strategy: { label: 'WGS', value: 'WGS' },
    signatureSetName: {
      label: 'COSMIC_v3_Signatures_GRCh37_SBS96',
      value: 'COSMIC_v3_Signatures_GRCh37_SBS96',
    },
    cancer: { label: 'Lung-AdenoCA', value: 'Lung-AdenoCA' },
  };

  const defaultValues = {
    ...publicDataDefaults,
    useAllCancer: false,
    compareEnabled: false,
    compare: publicDataDefaults,
  };

  const {
    control,
    handleSubmit,
    reset: resetForm,
    setValue,
    watch,
  } = useForm({ defaultValues });

  const { compareEnabled } = watch();

  function handleReset() {
    window.location.hash = '#/exploration';
    resetForm();
    dispatch(resetExplorationApi);
    resetExploration();
  }

  async function onSubmit({ compare, compareEnabled, ...data }) {
    try {
      mergeForm(data);
      mergeCompareForm(
        compareEnabled ? { ...compare, useAllCancer: data.useAllCancer } : {}
      );
      mergeMain({
        displayTab: 'tmb',
        submitted: true,
        loading: true,
        openSidebar: false,
        compare: compareEnabled,
      });
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
    mergeMain({ loading: false });
  }

  return (
    <Form onSubmit={handleSubmit(onSubmit)}>
      <LoadingOverlay active={isFetching} />
      {optionsError && <p>There was an error retrieving public data options</p>}
      <PublicDataFields
        control={control}
        setValue={setValue}
        watch={watch}
        options={options}
        disabled={submitted}
      />
      <hr className="mb-3" />
      <Form.Group controlId="useAllCancer" className="d-flex">
        <Form.Label className="mr-auto">
          Use all Cancer Types or Groups for TMB and MS Burden plots
        </Form.Label>
        <Form.Check inline id="useAllCancer">
          <Controller
            name="useAllCancer"
            control={control}
            render={({ field }) => (
              <Form.Check.Input
                {...field}
                type="checkbox"
                disabled={submitted}
                checked={field.value}
              />
            )}
          />
        </Form.Check>
      </Form.Group>
      <hr className="mb-3" />
      <Form.Group controlId="compareEnabled" className="d-flex">
        <Form.Label className="mr-auto">Public Data for Comparison</Form.Label>
        <Form.Check inline>
          <Controller
            name="compareEnabled"
            control={control}
            render={({ field }) => (
              <Form.Check.Input
                {...field}
                type="checkbox"
                disabled={submitted}
                checked={field.value}
              />
            )}
          />
        </Form.Check>
      </Form.Group>

      {compareEnabled && (
        <div className="border rounded p-2 mb-3 bg-light">
          <p className="font-weight-bold mb-2" style={{ fontSize: '12px' }}>
            Comparison Public Data
          </p>
          <PublicDataFields
            control={control}
            setValue={setValue}
            watch={watch}
            options={options}
            disabled={submitted}
            prefix="compare"
          />
        </div>
      )}

      <Row>
        <Col>
          <Button
            disabled={isFetching}
            className="w-100"
            variant="secondary"
            onClick={() => handleReset()}
          >
            Reset
          </Button>
        </Col>
        <Col>
          <Button
            disabled={isFetching || submitted}
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
