import Select from '@/components/controls/select/selectHookForm';
import { signatureSetCompare } from '@/services/utils';

// Study/Strategy/Signature Set/Cancer selects shared by the primary and compare Public Data segments.
export default function PublicDataFields({
  control,
  setValue,
  watch,
  options,
  disabled,
  prefix = '',
}) {
  const field = (name) => (prefix ? `${prefix}.${name}` : name);

  const { study, strategy, signatureSetName } = watch(
    prefix ? prefix : undefined
  ) || {};

  const studyOptions = options
    ? [...new Set(options.map((e) => e.study))].sort().map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  const strategyOptions = (study) => {
    if (options && study?.value) {
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
    if (options && study?.value && strategy?.value) {
      return [
        ...new Set(
          options
            .filter(
              (e) => e.study == study.value && e.strategy == strategy.value
            )
            .map((e) => e.signatureSetName)
        ),
      ]
        .sort(signatureSetCompare)
        .map((e) => ({ label: e, value: e }));
    } else {
      return [];
    }
  };

  const cancerOptions = (study, strategy, signatureSetName) => {
    if (
      options &&
      study?.value &&
      strategy?.value &&
      signatureSetName?.value
    ) {
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
    const signatureSets = signatureSetOptions(study, strategies[0]);
    const cancers = cancerOptions(study, strategies[0], signatureSets[0]);

    setValue(field('study'), study);
    setValue(field('strategy'), strategies[0]);
    setValue(field('signatureSetName'), signatureSets[0]);
    setValue(field('cancer'), cancers[0]);
  }

  function handleStrategyChange(strategy) {
    const signatureSets = signatureSetOptions(study, strategy);
    const cancers = cancerOptions(study, strategy, signatureSets[0]);

    setValue(field('strategy'), strategy);
    setValue(field('signatureSetName'), signatureSets[0]);
    setValue(field('cancer'), cancers[0]);
  }

  function handleSignatureSetChange(signatureSetName) {
    const cancers = cancerOptions(study, strategy, signatureSetName);

    setValue(field('signatureSetName'), signatureSetName);
    setValue(field('cancer'), cancers[0]);
  }

  return (
    <>
      <Select
        className="mb-2"
        name={field('study')}
        label="Study"
        disabled={disabled}
        options={studyOptions}
        control={control}
        onChange={handleStudyChange}
      />
      <Select
        className="mb-4"
        name={field('strategy')}
        label="Experimental Strategy"
        disabled={disabled}
        options={strategyOptions(study)}
        control={control}
        onChange={handleStrategyChange}
      />
      <Select
        className="mb-4"
        name={field('signatureSetName')}
        label="Reference Signature Set"
        disabled={disabled}
        options={signatureSetOptions(study, strategy)}
        control={control}
        onChange={handleSignatureSetChange}
      />
      <Select
        className="mb-2"
        name={field('cancer')}
        label="Cancer Type or Group"
        disabled={disabled}
        options={cancerOptions(study, strategy, signatureSetName)}
        control={control}
      />
    </>
  );
}
