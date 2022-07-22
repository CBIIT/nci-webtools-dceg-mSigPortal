import { useState } from 'react';
import { Form, Button, Row, Col } from 'react-bootstrap';
import { Controller, useForm } from 'react-hook-form';
import Select from '../../../controls/select/selectForm';
import { useSelector, useDispatch } from 'react-redux';
import { actions as exposureActions } from '../../../../services/store/exposure';
import { actions as modalActions } from '../../../../services/store/modal';
import {
  resetExplorationApi,
  useExplorationOptionsQuery,
} from '../../../../services/store/rootApi';
import {
  useExplorationSamplesMutation,
  useExplorationPublicMutation,
} from './apiSlice';

const actions = { ...exposureActions, ...modalActions };

export default function PublicForm() {
  const store = useSelector((state) => state.exposure);
  const { submitted } = store.main;

  const dispatch = useDispatch();
  const mergeForm = (state) =>
    dispatch(actions.mergeExposure({ publicForm: state }));
  const resetExposure = (_) => dispatch(actions.resetExposure());
  const mergeMain = (state) => dispatch(actions.mergeExposure({ main: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const {
    data: options,
    error: optionsError,
    isFetching: fetchingOptions,
  } = useExplorationOptionsQuery();

  const [fetchSamples, { isLoading, reset: resetSamples }] =
    useExplorationSamplesMutation();
  // const [handleSubmitWeb, { isLoading, reset: resetWeb }] =
  //   useExplorationPublicMutation();

  const defaultValues = {
    study: { label: 'PCAWG', value: 'PCAWG' },
    strategy: { label: 'WGS', value: 'WGS' },
    signatureSetName: {
      label: 'COSMIC_v3_Signatures_GRCh37_SBS96',
      value: 'COSMIC_v3_Signatures_GRCh37_SBS96',
    },
    cancer: { label: 'Lung-AdenoCA', value: 'Lung-AdenoCA' },
    cancerOnly: true,
  };

  const {
    control,
    handleSubmit,
    reset: resetForm,
    setValue,
    watch,
  } = useForm({ defaultValues });

  const { study, strategy, signatureSetName, cancer } = watch();

  function handleReset() {
    window.location.hash = '#/exploration';
    resetForm();
    // resetSamples();
    dispatch(resetExplorationApi);
    resetExposure();
  }

  async function onSubmit(data) {
    try {
      mergeMain({ submitted: true, loading: true });
      mergeForm(data);

      const params = {
        study: data.study.value,
        strategy: data.strategy.value,
        signatureSetName: data.signatureSetName.value,
      };
      const { samples, signatureNames } = await fetchSamples(params).unwrap();
      mergeMain({ displayTab: 'tmb', samples, signatureNames });
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

  const studyOptions = options
    ? [...new Set(options.map((e) => e.study))].map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  const strategyOptions = (study) => {
    if (options && study.value) {
      return [
        ...new Set(
          options.filter((e) => e.study == study.value).map((e) => e.strategy)
        ),
      ].map((e) => ({
        label: e,
        value: e,
      }));
    } else {
      return [];
    }
  };

  const signatureSetOptions = (study, strategy) => {
    if (options && study.value && strategy.value) {
      return [
        ...new Set(
          options
            .filter(
              (e) => e.study == study.value && e.strategy == strategy.value
            )
            .map((e) => e.signatureSetName)
        ),
      ].map((e) => ({
        label: e,
        value: e,
      }));
    } else {
      return [];
    }
  };

  // const getSignatureNameOptions = () => {
  //   if (results) {
  //     return [...new Set(results.map((e) => e.signatureName))].map((e) => ({
  //       label: e,
  //       value: e,
  //     }));
  //   } else {
  //     return [];
  //   }
  // };

  // const signatureNameOptions = getSignatureNameOptions();

  const cancerOptions = (study, strategy, signatureSetName) => {
    if (options && study.value && strategy.value && signatureSetName.value) {
      return [
        ...new Set(
          options
            .filter(
              (e) =>
                e.study == study.value &&
                e.strategy == strategy.value &&
                e.signatureSetName == signatureSetName.value
            )
            .map((e) => e.cancer)
        ),
      ].map((e) => ({
        label: e,
        value: e,
      }));
    } else {
      return [];
    }
  };

  function handleStudyChange(study) {
    const strategies = strategyOptions(study);
    const signatureSets = signatureSetOptions(study, strategies[0]);
    const cancers = cancerOptions(study, strategies[0], signatureSets[0]);

    setValue('study', study);
    setValue('strategy', strategies[0]);
    setValue('signatureSetName', signatureSets[0]);
    setValue('cancer', cancers[0]);
  }

  function handleStrategyChange(strategy) {
    const signatureSets = signatureSetOptions(study, strategy);
    const cancers = cancerOptions(study, strategy, signatureSets[0]);

    setValue('strategy', strategy);
    setValue('signatureSetName', signatureSets[0]);
    setValue('cancer', cancers[0]);
  }

  function handleSignatureSetChange(signatureSetName) {
    const cancers = cancerOptions(study, strategy, signatureSetName);

    setValue('signatureSetName', signatureSetName);
    setValue('cancer', cancers[0]);
  }

  return (
    <Form onSubmit={handleSubmit(onSubmit)}>
      {optionsError && <p>There was an error retrieving public data options</p>}
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
        className="mb-4"
        name="strategy"
        label="Experimental Strategy"
        disabled={submitted}
        options={strategyOptions(study)}
        control={control}
        onChange={handleStrategyChange}
      />
      <Select
        className="mb-4"
        name="signatureSetName"
        label="Reference Signature Set"
        disabled={submitted}
        options={signatureSetOptions(study, strategy)}
        control={control}
        onChange={handleSignatureSetChange}
      />
      <Select
        className="mb-2"
        name="cancer"
        label="Cancer Type or Group"
        disabled={submitted}
        options={cancerOptions(study, strategy, signatureSetName)}
        control={control}
      />
      <Form.Group>
        <Controller
          name="cancerOnly"
          control={control}
          render={({ field }) => (
            <Form.Check
              {...field}
              id="cancerOnly"
              disabled={submitted}
              type="checkbox"
              label="Cancer Type or Group Only"
              checked={field.value}
            />
          )}
        />
      </Form.Group>

      <Row>
        <Col>
          <Button
            disabled={fetchingOptions || isLoading}
            className="w-100"
            variant="secondary"
            onClick={() => handleReset()}
          >
            Reset
          </Button>
        </Col>
        <Col>
          <Button
            disabled={fetchingOptions || isLoading || submitted}
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
