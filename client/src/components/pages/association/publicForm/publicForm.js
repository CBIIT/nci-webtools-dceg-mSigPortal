import { useEffect, useState } from 'react';
import { useForm } from 'react-hook-form';
import { Form, Row, Col, Button } from 'react-bootstrap';
import axios from 'axios';
import Select from '../../../controls/select/selectForm';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { useSelector, useDispatch } from 'react-redux';
import { actions as associationActions } from '../../../../services/store/association';
import { actions as modalActions } from '../../../../services/store/modal';
import {
  resetAssociationApi,
  useAssociationOptionsQuery,
  useExposureOptionsQuery,
} from '../../../../services/store/rootApi';
import { useAssociationPublicMutation } from './apiSlice';

const actions = { ...associationActions, ...modalActions };

export default function PublicForm() {
  const { submitted } = useSelector((state) => state.association.main);

  const dispatch = useDispatch();
  const mergeState = async (state) =>
    dispatch(actions.mergeAssociation({ main: state }));
  const resetAssociation = (_) => dispatch(actions.resetAssociation());
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const defaultValues = {
    study: { label: 'PCAWG', value: 'PCAWG' },
    strategy: { label: 'WGS', value: 'WGS' },
    rsSet: {
      label: 'COSMIC_v3_Signatures_GRCh37_SBS96',
      value: 'COSMIC_v3_Signatures_GRCh37_SBS96',
    },
    cancer: { label: 'Lung-AdenoCA', value: 'Lung-AdenoCA' },
  };

  const {
    control,
    handleSubmit,
    reset: resetForm,
    setValue,
    watch,
  } = useForm({ defaultValues });

  const formStudy = watch('study');
  const formStrategy = watch('strategy');
  const formSignatureSet = watch('rsSet');

  const [loading, setLoading] = useState(false);

  const {
    data: associationOptions,
    error: associationError,
    isFetching: fetchingAssociation,
  } = useAssociationOptionsQuery();
  const {
    data: exposureOptions,
    error: exposureError,
    isFetching: fetchingExposure,
  } = useExposureOptionsQuery({
    study: formStudy?.value,
    strategy: formStrategy?.value,
  });
  const [submitAssociation, { isLoading }] = useAssociationPublicMutation();

  function handleReset() {
    resetForm();
    dispatch(resetAssociationApi);
    resetAssociation();
  }

  async function onSubmit(data) {
    console.log(data);
    const params = {
      study: data.study.value,
      strategy: data.strategy.value,
      rsSet: data.rsSet.value,
      cancer: data.cancer.value,
    };
    const getAssocVarData = () =>
      axios
        .post('web/associationWrapper', {
          fn: 'getAssocVarData',
          args: params,
        })
        .then(({ data }) => data);

    const getExpVarData = () =>
      axios
        .post('web/associationWrapper', {
          fn: 'getExpVarData',
          args: params,
        })
        .then(({ data }) => data);

    try {
      setLoading(true);
      mergeState({ submitted: true, ...data });

      const [assocResponse, expResponse] = await Promise.all([
        getAssocVarData(),
        getExpVarData(),
      ]);

      const { id, output: assocOutput } = assocResponse;
      const { output: expOutput } = expResponse;

      if (assocOutput.uncaughtError) throw assocOutput.uncaughtError;
      if (expOutput.uncaughtError) throw expOutput.uncaughtError;

      mergeState({
        submitted: true,
        assocVarData: assocOutput.assocVarData,
        assocFullDataPath: assocOutput.fullDataPath,
        expVarList: expOutput.expVarList,
        displayTab: 'univariable',
        openSidebar: false,
      });
      dispatch(actions.mergeAssociation({ univariable: { id } }));
    } catch (error) {
      if (error.originalStatus == 504) {
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
        mergeError(error.data);
      }
    } finally {
      setLoading(false);
    }
  }

  // automatically select new signature set and cancer when new exposure options are retrieved after study changes
  useEffect(() => {
    if (exposureOptions) {
      const signatureSets = signatureSetOptions(formStudy, formStrategy);
      handleSignatureSetChange(signatureSets[0]);
    }
  }, [exposureOptions]);

  const studyOptions = associationOptions
    ? [...new Set(associationOptions.map((e) => e.study))]
        .sort(
          new Intl.Collator('en', { numeric: true, sensitivity: 'accent' })
            .compare
        )
        .map((e) => ({
          label: e,
          value: e,
        }))
    : [];

  const strategyOptions = (study) => {
    if (associationOptions && study.value) {
      return [
        ...new Set(
          associationOptions
            .filter((e) => e.study == study.value)
            .map((e) => e.strategy)
        ),
      ]
        .sort()
        .map((e) => ({
          label: e,
          value: e,
        }));
    } else {
      return [];
    }
  };

  const signatureSetOptions = (study, strategy) => {
    if (exposureOptions && study.value && strategy.value) {
      return [
        ...new Set(
          exposureOptions
            .filter(
              (e) => e.study == study.value && e.strategy == strategy.value
            )
            .map((e) => e.signatureSetName)
        ),
      ]
        .sort()
        .map((e) => ({
          label: e,
          value: e,
        }));
    } else {
      return [];
    }
  };

  const cancerOptions = (study, strategy, signatureSet) => {
    if (
      exposureOptions &&
      study.value &&
      strategy.value &&
      signatureSet.value
    ) {
      return [
        ...new Set(
          exposureOptions
            .filter(
              (e) =>
                e.study == study.value &&
                e.strategy == strategy.value &&
                e.signatureSetName == signatureSet.value
            )
            .sort()
            .map((e) => e.cancer)
        ),
      ]
        .sort()
        .map((e) => ({
          label: e,
          value: e,
        }));
    } else {
      return [];
    }
  };

  function handleStudyChange(study) {
    const strategies = strategyOptions(study);

    setValue('study', study);
    setValue('strategy', strategies[0]);
  }

  function handleSignatureSetChange(signatureSet) {
    const cancers = cancerOptions(formStudy, formStrategy, signatureSet);

    setValue('rsSet', signatureSet);
    setValue('cancer', cancers[0]);
  }

  return (
    <Form onSubmit={handleSubmit(onSubmit)}>
      <LoadingOverlay
        active={fetchingAssociation || fetchingExposure || loading}
      />
      {(associationError || exposureError) && (
        <p>There was an error retrieving public data options</p>
      )}
      <Select
        className="mb-2"
        name="study"
        label="Study"
        disabled={submitted || fetchingAssociation}
        options={studyOptions}
        control={control}
        onChange={handleStudyChange}
      />
      <Select
        className="mb-2"
        name="strategy"
        label="Experimental Strategy"
        disabled={submitted || fetchingAssociation}
        options={strategyOptions(formStudy)}
        control={control}
      />
      <Select
        className="mb-2"
        name="rsSet"
        label="Reference Signature Set"
        disabled={submitted || fetchingExposure}
        options={signatureSetOptions(formStudy, formStrategy)}
        control={control}
        onChange={handleSignatureSetChange}
      />
      <Select
        className="mb-4"
        name="cancer"
        label="Cancer Type or Group"
        disabled={submitted || fetchingAssociation}
        options={cancerOptions(formStudy, formStrategy, formSignatureSet)}
        control={control}
      />
      <Row>
        <Col>
          <Button
            disabled={fetchingAssociation || fetchingExposure}
            className="w-100"
            variant="secondary"
            onClick={() => handleReset()}
          >
            Reset
          </Button>
        </Col>
        <Col>
          <Button
            disabled={fetchingAssociation || fetchingExposure || submitted}
            className="w-100"
            variant="primary"
            type="submit"
          >
            Load Study
          </Button>
        </Col>
      </Row>
    </Form>
  );
}
