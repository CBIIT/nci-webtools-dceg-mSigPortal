import { useState, useEffect } from 'react';
import { Form, Button, Row, Col } from 'react-bootstrap';
import { useForm, Controller } from 'react-hook-form';
import { useHistory, useParams } from 'react-router-dom';
import SelectForm from '../../controls/select/selectHookForm';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { useSelector, useDispatch } from 'react-redux';
import { actions as extractionActions } from '../../../services/store/extraction';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faGithub } from '@fortawesome/free-brands-svg-icons';
import { actions as modalActions } from '../../../services/store/modal';
import {
  resetExtractionApi,
  useSeqmatrixOptionsQuery,
  useSignatureOptionsQuery,
  useRefGenomeQuery,
} from '../../../services/store/rootApi';
import {
  useUploadMutation,
  useSubmitMutation,
  useParamsQuery,
} from './apiSlice';

const actions = { ...extractionActions, ...modalActions };

export default function ExtractionForm({ formLimits }) {
  const { submitted, ...state } = useSelector((state) => state.extraction);
  const id = useParams().id || state.id || false;
  const history = useHistory();

  const dispatch = useDispatch();
  const resetExtraction = (_) => dispatch(actions.resetExtraction());
  const mergeState = (state) => dispatch(actions.mergeExtraction(state));
  const mergeSuccess = (msg) =>
    dispatch(actions.mergeModal({ success: { visible: true, message: msg } }));

  const { data: params } = useParamsQuery(id, { skip: !id });

  const [selectedOptions, setSelectedOptions] = useState([]);
  const [isDefaultContext, setIsDefaultContext] = useState(true); // Set the default value according to your requirements

  const [warning, setWarning] = useState('');

  // query options to populate form
  const {
    data: seqmatrixOptions,
    error: seqmatrixError,
    isFetching: fetchingSeqmatrixOptions,
  } = useSeqmatrixOptionsQuery({
    columns: ['study', 'strategy', 'cancer', 'profile', 'matrix'],
  });
  const {
    data: signatureOptions,
    error: signatureError,
    isFetching: fetchingSignatureOptions,
  } = useSignatureOptionsQuery();
  const {
    data: genomeData,
    error: genomeError,
    isFetching: fetchingGenomes,
  } = useRefGenomeQuery();

  const [uploadFiles, { isLoading: loadingUpload }] = useUploadMutation();
  const [submitForm, { isLoading: loadingSubmit }] = useSubmitMutation();

  // toggle visibility of advanced menu
  const [showAdvanced, setShowAdvanced] = useState(false);

  // define various select options
  const genomeOptions = genomeData
    ? [...new Set(genomeData.map((e) => e.genome))].sort().map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  const signatureNameOptions = (signatureSetName) =>
    signatureOptions && signatureSetName
      ? [
          { label: 'all', value: 'all' },
          ...[
            ...new Set(
              signatureOptions
                .filter((e) => e.signatureSetName == signatureSetName.value)
                .map((e) => e.signatureName)
                .sort((a, b) =>
                  a.localeCompare(b, undefined, {
                    numeric: true,
                    sensitivity: 'base',
                  })
                )
            ),
          ].map((e) => ({ label: e, value: e })),
        ]
      : [];

  const [filteredSignatureSetOptions, setFilteredSignatureSetOptions] =
    useState();

  const signatureSetOptions = signatureOptions
    ? [...new Set(signatureOptions.map((e) => e.signatureSetName))]
        .map((e) => ({ label: e, value: e }))
        .sort((a, b) => a.value.localeCompare(b.value))
    : [];

  const handleContextTypeChange = (selectedOption) => {
    setValue('context_type', selectedOption);
    let filteredOptions;
    if (selectedOption.value !== 'default') {
      setIsDefaultContext(false);
      filteredOptions = signatureSetOptions
        .filter((option) => option.value.includes(selectedOption.value))
        .sort((a, b) => a.value.localeCompare(b.value));
    } else {
      setIsDefaultContext(true);
      filteredOptions = signatureSetOptions;
    }
    setFilteredSignatureSetOptions(filteredOptions);

    if (selectedOption.value === 'SBS96') {
      const cosmicOption = filteredOptions.find(
        (option) => option.value === 'COSMIC_v3.3_Signatures_GRCh37_SBS96'
      );
      setValue('signatureSetName', cosmicOption);
    } else {
      setValue('signatureSetName', filteredOptions[0]);
    }
  };
  const dataTypeOptions = ['matrix'].map((e) => ({
    label: e,
    value: e,
  }));

  const studyOptions = seqmatrixOptions
    ? [...new Set(seqmatrixOptions.map((e) => e.study))].sort().map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  const cancerOptions = (study) => {
    if (seqmatrixOptions && study.value) {
      const options = [
        ...[
          ...new Set(
            seqmatrixOptions
              .filter((e) => e.study == study.value)
              .map((e) => e.cancer)
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
    if (seqmatrixOptions && study.value) {
      const options =
        cancer.value == '*ALL'
          ? [
              ...new Set(
                seqmatrixOptions
                  .filter((e) => e.study == study.value)
                  .map((e) => e.strategy)
              ),
            ]
          : [
              ...new Set(
                seqmatrixOptions
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

  // define onChange handlers
  function handleStudyChange(study) {
    const cancers = cancerOptions(study);
    const strategies = strategyOptions(study, cancers[0]);

    setValue('study', study);
    setValue('cancer', cancers[0]);
    setValue('strategy', strategies[0]);
  }

  function handleCancerChange(cancer) {
    const strategies = strategyOptions(study, cancer);

    setValue('cancer', cancer);
    setValue('strategy', strategies[0]);
  }

  // define form
  const defaultValues = {
    source: 'public',
    study: { label: 'PCAWG', value: 'PCAWG' },
    cancer: { label: 'Lung-AdenoCA', value: 'Lung-AdenoCA' },
    strategy: { label: 'WGS', value: 'WGS' },
    input_type: { label: 'matrix', value: 'matrix' },
    reference_genome: { label: 'GRCh37', value: 'GRCh37' },
    opportunity_genome: { label: 'GRCh37', value: 'GRCh37' },
    exome: false,
    signatureSetName: {
      label: 'COSMIC_v3.3_Signatures_GRCh37_SBS96',
      value: 'COSMIC_v3.3_Signatures_GRCh37_SBS96',
    },
    signatureName: [{ label: 'all', value: 'all' }],
    extractTool: {
      label: 'SigProfilerExtractor',
      value: 'SigProfilerExtractor',
    },

    seeds: 'random',
    matrix_normalization: 'gmm',
    nmf_init: 'random',
    precision: 'single',

    minimum_signatures: 1,
    maximum_signatures: 12,
    nmf_replicates: 50,
    min_nmf_iterations: 1000,
    max_nmf_iterations: 100000,
    nmf_test_conv: 10000,

    email: '',
    jobName: '',
  };

  const sbsSample = {
    ...defaultValues,
    source: 'user',
    context_type: { label: 'SBS96', value: 'SBS96' },
    minimum_signatures: 1,
    maximum_signatures: 2,
    min_nmf_iterations: 2,
    max_nmf_iterations: 4,
    nmf_test_conv: 2,
    nmf_replicates: 10,
  };
  const dbsSample = {
    ...defaultValues,
    source: 'user',
    context_type: { label: 'DBS78', value: 'DBS78' },
    signatureSetName: {
      label: 'COSMIC_v3.3_Signatures_GRCh37_DBS78',
      value: 'COSMIC_v3.3_Signatures_GRCh37_DBS78',
    },
    minimum_signatures: 1,
    maximum_signatures: 2,
    min_nmf_iterations: 2,
    max_nmf_iterations: 4,
    nmf_test_conv: 2,
    nmf_replicates: 10,
  };

  const {
    control,
    register,
    handleSubmit,
    reset: resetForm,
    setValue,
    watch,
    formState: { errors },
  } = useForm({ defaultValues: defaultValues });

  const {
    source,
    study,
    cancer,
    inputFile,
    reference_genome,
    context_type,
    signatureSetName,
    input_type,
  } = watch();

  const contextTypeOptions = (() => {
    if (source === 'user') {
      return signatureOptions
        ? [
            { label: 'default', value: 'default' },
            ...[...new Set(signatureOptions.map((e) => e.profile + e.matrix))]
              .sort((a, b) => a.localeCompare(b, 'en', { numeric: true }))
              .map((e) => ({ label: e, value: e })),
          ]
        : [];
    } else {
      if (signatureOptions && seqmatrixOptions) {
        const sigProfiles = [
          ...new Set(signatureOptions.map((e) => e.profile + e.matrix)),
        ];
        const seqmatrixProfiles = [
          ...new Set(seqmatrixOptions.map((e) => e.profile + e.matrix)),
        ];
        const commonProfiles = sigProfiles.filter((e) =>
          seqmatrixProfiles.includes(e)
        );
        return [
          { label: 'default', value: 'default' },
          ...commonProfiles
            .sort((a, b) => a.localeCompare(b, 'en', { numeric: true }))
            .map((e) => ({ label: e, value: e })),
        ];
      } else {
        return [];
      }
    }
  })();

  // update url with id
  // useEffect(() => {
  //   if (!id && state.id) history.push(`/extraction/${state.id}`);
  // }, [id, state.id]);

  // populate form
  useEffect(() => {
    if (params && id) {
      resetForm(params.form);
      mergeState({ submitted: true, id });
    }
  }, [params, id]);

  // set inital genome option after querying options
  useEffect(() => {
    if (!reference_genome && genomeOptions.length)
      setValue('reference_genome', genomeOptions[0]);
  }, [genomeOptions]);
  // set inital context type option after querying options
  useEffect(() => {
    if (!context_type && contextTypeOptions.length)
      setValue('context_type', contextTypeOptions[0]);
  }, [contextTypeOptions]);

  function handleReset() {
    history.push('/extraction');
    resetForm(defaultValues);
    resetExtraction();
    dispatch(resetExtractionApi);
  }

  async function onSubmit(data) {
    mergeState({ submitted: true });
    const formData = new FormData();
    formData.append('inputFile', data.inputFile);
    const { id } = await uploadFiles(formData).unwrap();

    const args = {
      ...(source === 'user' && {
        input_type: data.input_type.value,
        input_data: data.inputFile.name,
      }),
      ...(source === 'public' && {
        input_type: 'matrix',
        //input_data: 'ExtractionData.all',
      }),

      reference_genome: data.reference_genome.value,
      opportunity_genome: data.opportunity_genome.value,
      exome: data.exome ? 'True' : 'False',
      context_type: data.context_type.value,
      minimum_signatures: data.minimum_signatures,
      maximum_signatures: data.maximum_signatures,
      nmf_replicates: data.nmf_replicates,
      resample: data.resample ? 'True' : 'False',
      seeds: data.seeds,
      min_nmf_iterations: data.min_nmf_iterations,
      max_nmf_iterations: data.max_nmf_iterations,
      nmf_test_conv: data.nmf_test_conv,
      stability: data.stability,
      min_stability: data.min_stability,
      combined_stability: data.combined_stability,
      allow_stability_drop: data.allow_stability_drop ? 'True' : 'False',
    };

    const signatureSetName = data.signatureSetName.value;
    const profile =
      data.context_type.value == 'default'
        ? signatureOptions.filter(
            (e) => e.signatureSetName == signatureSetName
          )[0].profile
        : data.context_type.value.match(/^\D*/)[0];
    const matrix =
      data.context_type.value == 'default'
        ? signatureOptions.filter(
            (e) => e.signatureSetName == signatureSetName
          )[0].matrix + ''
        : data.context_type.value.match(/\d*$/)[0];
    const signatureQuery = {
      signatureSetName,
      profile,
      matrix,
      ...(data.signatureName[0].value != 'all' && {
        signatureName: data.signatureName.map((e) => e.value).join(';'),
      }),
    };

    const seqmatrixQuery = {
      study: data.study?.value,
      cancer: data.cancer?.value,
      strategy: data.strategy?.value,
      profile,
      matrix,
    };

    const params = {
      args,
      signatureQuery,
      seqmatrixQuery,
      id,
      email: data.email,
      jobName: data.jobName || 'Extraction',
      form: data,
    };
    //console.log('Params----', params);
    const submitStatus = await submitForm(params).unwrap();

    history.push(`/extraction/${submitStatus.id}`);
    mergeState({ id, displayTab: 'status' });
    mergeSuccess(
      `Most Jobs take a long time, you will receive an email when the extraction job is complete. It is safe to close the window now`
    );
    const jobs = JSON.parse(localStorage.getItem('jobs')) || [];
    localStorage.setItem('jobs', JSON.stringify([...jobs, id]));
  }

  const handleSelectChange = (selectedOptions) => {
    const isAllOptionSelected = selectedOptions.some(
      (option) => option.value === 'all'
    );

    if (isAllOptionSelected) {
      // Clear the warning when 'all' is selected and unselect individual signatures
      setWarning('');
      selectedOptions = [{ value: 'all', label: 'all' }];
    } else {
      // Filter out the 'all' option if it exists in selectedOptions
      selectedOptions = selectedOptions.filter(
        (option) => option.value !== 'all'
      );

      if (selectedOptions.length <= 1) {
        // Show a warning message if only one option (other than 'all') is selected
        setWarning('More than one signature is required as input');
      } else {
        // Clear the warning if multiple options are selected
        setWarning('');
      }
    }

    // Update the form values with the modified selectedOptions
    setValue('signatureName', selectedOptions);
  };

  return (
    <div className="p-3 bg-white border rounded">
      <Form onSubmit={handleSubmit(onSubmit)}>
        <LoadingOverlay
          active={fetchingSeqmatrixOptions || loadingUpload || loadingSubmit}
        />
        {(seqmatrixError || signatureError) && (
          <p>There was an error retrieving public data options</p>
        )}
        <div style={{ maxHeight: '900px', overflow: 'hidden auto' }}>
          <div className="border rounded p-2 mb-3">
            <Form.Group>
              <Form.Label className="mr-4">Data Source</Form.Label>
              <Controller
                name="source"
                control={control}
                render={({ field }) => (
                  <Form.Check
                    {...field}
                    inline
                    id="radioPublic"
                    type="radio"
                    label={<span className="font-weight-normal">Public</span>}
                    value={'public'}
                    checked={field.value === 'public'}
                    disabled={submitted || id}
                  />
                )}
              />
              <Controller
                name="source"
                control={control}
                render={({ field }) => (
                  <Form.Check
                    {...field}
                    inline
                    id="radioUser"
                    type="radio"
                    label={<span className="font-weight-normal">User</span>}
                    value={'user'}
                    checked={field.value === 'user'}
                    disabled={submitted || id}
                  />
                )}
              />
            </Form.Group>
            {source === 'public' ? (
              <div key="public">
                <SelectForm
                  className="mb-2"
                  name="study"
                  label="Study"
                  disabled={submitted || id || fetchingSeqmatrixOptions}
                  options={studyOptions}
                  control={control}
                  onChange={handleStudyChange}
                />
                <SelectForm
                  className="mb-2"
                  name="cancer"
                  label="Cancer Type or Group"
                  disabled={submitted || id || fetchingSeqmatrixOptions}
                  options={cancerOptions(study)}
                  control={control}
                  onChange={handleCancerChange}
                />
                <SelectForm
                  className="mb-2"
                  name="strategy"
                  label="Experimental Strategy"
                  disabled={submitted || id || fetchingSeqmatrixOptions}
                  options={strategyOptions(study, cancer)}
                  control={control}
                />
              </div>
            ) : (
              <div key="user">
                <SelectForm
                  name="input_type"
                  label="Data Type"
                  disabled={submitted || id}
                  options={dataTypeOptions}
                  control={control}
                />
                <Form.Group>
                  <Form.Label>
                    Upload File <span style={{ color: 'crimson' }}>*</span>
                  </Form.Label>
                  <Controller
                    name="inputFile"
                    control={control}
                    rules={{ required: source == 'user' }}
                    render={({ field }) => (
                      <Form.File
                        {...field}
                        value={''} // set dummy value for file input
                        disabled={submitted || id}
                        id="inputFile"
                        label={
                          inputFile?.name ||
                          params?.args.input_data ||
                          'Upload Data File...'
                        }
                        isInvalid={errors.inputFile}
                        feedback="Please upload a data file"
                        onChange={(e) => {
                          if (e.target.files.length) {
                            setValue('inputFile', e.target.files[0]);
                          }
                        }}
                        custom
                      />
                    )}
                  />
                </Form.Group>
                <Button
                  variant="link"
                  disabled={submitted || id}
                  onClick={async () => {
                    resetForm(sbsSample);
                    const file = 'extraction_sample_SBS96.all';
                    const path = 'assets/exampleInput/' + file;
                    setValue(
                      'inputFile',
                      new File([await (await fetch(path)).blob()], file)
                    );
                  }}
                >
                  Load SBS Matrix
                </Button>
                <Button
                  variant="link"
                  disabled={submitted || id}
                  onClick={async () => {
                    resetForm(dbsSample);
                    const file = 'extraction_sample_DBS78.all';
                    const path = 'assets/exampleInput/' + file;
                    setValue(
                      'inputFile',
                      new File([await (await fetch(path)).blob()], file)
                    );
                  }}
                >
                  Load DBS Matrix
                </Button>
              </div>
            )}
          </div>

          <SelectForm
            name="reference_genome"
            label="Reference Genome Build"
            disabled={submitted || id}
            options={genomeOptions}
            control={control}
            rules={{ required: true }}
          />
          <SelectForm
            name="opportunity_genome"
            label="Opportunity Genome Build"
            disabled={submitted || id}
            options={genomeOptions}
            control={control}
          />
          <Form.Group>
            <Controller
              name="exome"
              control={control}
              render={({ field }) => (
                <Form.Check
                  {...field}
                  id="exome"
                  type="checkbox"
                  label={'Exome'}
                  checked={field.checked}
                  disabled={submitted || id}
                />
              )}
            />
          </Form.Group>
          <SelectForm
            name="context_type"
            label="Context Type"
            disabled={submitted || id}
            options={contextTypeOptions}
            control={control}
            onChange={handleContextTypeChange}
            rules={{ required: true }}
          />
          <SelectForm
            name="signatureSetName"
            label="Reference Signature Set"
            disabled={submitted || id}
            //options={signatureSetOptions}
            //options={filteredSignatureSetOptions}
            options={
              isDefaultContext
                ? signatureSetOptions
                : filteredSignatureSetOptions
            }
            control={control}
          />
          {warning && <div style={{ color: 'red' }}>{warning}</div>}
          <SelectForm
            name="signatureName"
            label="Select Signature Names"
            disabled={submitted || id}
            options={signatureNameOptions(signatureSetName)}
            control={control}
            onChange={handleSelectChange}
            // onChange={(selectedOptions) => {
            //   // Check if the 'all' option is selected
            //   const isAllOptionSelected = selectedOptions.some(
            //     (option) => option.value === 'all'
            //   );
            //   if (!isAllOptionSelected) {
            //     // Filter out the 'all' option if it exists in selectedOptions
            //     selectedOptions = selectedOptions.filter(
            //       (option) => option.value !== 'all'
            //     );

            //     if (selectedOptions.length === 1) {
            //       // Show a warning message if only one option (other than 'all') is selected
            //       setWarning('More than one signature is required as input');
            //     } else {
            //       // Clear the warning if multiple options are selected
            //       setWarning('');
            //     }
            //   } else {
            //     // Clear the warning when 'all' is selected
            //     setWarning('');
            //     setValue('signatureName', [{ value: 'all', label: 'all' }]);
            //   }

            //   // Update the form values with the modified selectedOptions
            //   setValue('signatureName', selectedOptions);
            // }}
            // onChange={(values, e) => {
            //   // remove "all" option if a specific signature is selected
            //   // remove other options if "all" is selected
            //   console.log('e', e);
            //   if (e.action === 'select-option') {
            //     if (e.option.value !== 'all') {
            //       setValue(
            //         'signatureName',
            //         values.filter((e) => e.value !== 'all')
            //       );
            //     } else if (e.option.value === 'all') {
            //       setValue('signatureName', [e.option]);
            //     }
            //   } else setValue('signatureName', values);
            // }}
            isMulti
          />
          <SelectForm
            name="extractTool"
            label={
              <div>
                Extract Tool{' '}
                <Button
                  variant="outline-secondary"
                  size="sm"
                  href="https://github.com/AlexandrovLab/SigProfilerExtractor"
                  target="_blank"
                  rel="noopener noreferrer"
                  style={{ color: 'black' }}
                  className="ml-2 mb-2"
                >
                  <FontAwesomeIcon icon={faGithub} className="mr-1" />
                  GitHub
                </Button>
              </div>
            }
            disabled={submitted || id}
            options={[
              { label: 'SigProfilerExtractor', value: 'SigProfilerExtractor' },
            ]}
            control={control}
          />

          <Button
            className="p-0"
            variant="link"
            onClick={(_) => setShowAdvanced(!showAdvanced)}
          >
            Advanced Parameters {showAdvanced ? '-' : '+'}
          </Button>

          <div className={showAdvanced ? 'd-block' : 'd-none'}>
            <fieldset className="border rounded p-2 mb-3">
              <legend className="font-weight-bold">NMF Replicates</legend>
              <Form.Group controlId="minSignatures">
                <Form.Label>Minimum Signatures</Form.Label>
                <Form.Control
                  {...register('minimum_signatures')}
                  type="number"
                  onWheel={(e) => e.target.blur()}
                  min={formLimits.minimum_signatures[0]}
                  max={formLimits.maximum_signatures[1]}
                  defaultValue={1}
                  disabled={submitted || id}
                />
              </Form.Group>
              <Form.Group controlId="maxSignatures">
                <Form.Label>Maximum Signatures</Form.Label>
                <Form.Control
                  {...register('maximum_signatures')}
                  type="number"
                  onWheel={(e) => e.target.blur()}
                  min={formLimits.maximum_signatures[0]}
                  max={formLimits.maximum_signatures[1]}
                  defaultValue={12}
                  disabled={submitted || id}
                />
              </Form.Group>
              <Form.Group controlId="nmf_replicates">
                <Form.Label>NMF Replicates Size</Form.Label>
                <Form.Control
                  {...register('nmf_replicates')}
                  type="number"
                  onWheel={(e) => e.target.blur()}
                  min={formLimits.nmf_replicates[0]}
                  max={formLimits.nmf_replicates[1]}
                  defaultValue={50}
                  disabled={submitted || id}
                />
              </Form.Group>
              <Form.Group controlId="resample">
                <Form.Check
                  {...register('resample')}
                  label="Resample"
                  id="resample"
                  defaultChecked={true}
                  disabled={submitted || id}
                />
              </Form.Group>
              <Form.Group controlId="seeds">
                <Form.Label>Seeds</Form.Label>
                <Form.Control
                  {...register('seeds')}
                  defaultValue="random"
                  disabled={submitted || id}
                />
              </Form.Group>
            </fieldset>

            <fieldset className="border rounded p-2 mb-3">
              <legend className="font-weight-bold">NMF Engines</legend>
              <Form.Group controlId="matrixNormalization">
                <Form.Label>Matrix Normalization</Form.Label>
                <Controller
                  name="matrix_normalization"
                  control={control}
                  render={({ field }) => (
                    <Form.Control
                      {...field}
                      as="select"
                      disabled={submitted || id}
                    >
                      <option>gmm</option>
                      <option>log2</option>
                      <option>custom</option>
                      <option>none</option>
                      {['gmm', 'log2', 'custom', 'none'].map((e) => (
                        <option value={e}>{e}</option>
                      ))}
                    </Form.Control>
                  )}
                />
              </Form.Group>
              <Form.Group controlId="nmfInit">
                <Form.Label>NMF Initialization</Form.Label>
                <Controller
                  name="nmf_init"
                  control={control}
                  render={({ field }) => (
                    <Form.Control
                      {...field}
                      as="select"
                      disabled={submitted || id}
                    >
                      {[
                        'random',
                        'nndsvd',
                        'nndsvda',
                        'nndsvdar',
                        'nndsvd_min',
                      ].map((e) => (
                        <option value={e}>{e}</option>
                      ))}
                    </Form.Control>
                  )}
                />
              </Form.Group>
              <Form.Group controlId="precision">
                <Form.Label>Precision</Form.Label>
                <Controller
                  name="precision"
                  control={control}
                  render={({ field }) => (
                    <Form.Control
                      {...field}
                      as="select"
                      disabled={submitted || id}
                    >
                      {['single', 'double'].map((e) => (
                        <option value={e}>{e}</option>
                      ))}
                    </Form.Control>
                  )}
                />
              </Form.Group>
              <Form.Group controlId="minNmfIterations">
                <Form.Label>Minimum NMF Iterations</Form.Label>
                <Form.Control
                  {...register('min_nmf_iterations', {})}
                  type="number"
                  onWheel={(e) => e.target.blur()}
                  min={formLimits.min_nmf_iterations[0]}
                  max={formLimits.min_nmf_iterations[1]}
                  defaultValue={1000}
                  disabled={submitted || id}
                />
              </Form.Group>
              <Form.Group controlId="maxNmfIterations">
                <Form.Label>Maximum NMF Iterations</Form.Label>
                <Form.Control
                  {...register('max_nmf_iterations')}
                  type="number"
                  onWheel={(e) => e.target.blur()}
                  min={formLimits.max_nmf_iterations[0]}
                  max={formLimits.max_nmf_iterations[1]}
                  defaultValue={100000}
                  disabled={submitted || id}
                />
              </Form.Group>
              <Form.Group controlId="nmfTestConv">
                <Form.Label>NMF Test Convergence</Form.Label>
                <Form.Control
                  {...register('nmf_test_conv')}
                  type="number"
                  onWheel={(e) => e.target.blur()}
                  defaultValue={1000}
                  min={formLimits.nmf_test_conv[0]}
                  max={formLimits.nmf_test_conv[1]}
                  disabled={submitted || id}
                />
              </Form.Group>
              <Form.Group controlId="nmfTolerance">
                <Form.Label>NMF Tolerance</Form.Label>
                <Form.Control
                  {...register('nmf_tolerance')}
                  type="number"
                  onWheel={(e) => e.target.blur()}
                  defaultValue={1.0e-15}
                  disabled={submitted || id}
                />
              </Form.Group>
            </fieldset>

            <fieldset className="border rounded p-2 mb-3">
              <legend className="font-weight-bold">
                Solution Estimation Thresholds
              </legend>
              <Form.Group>
                <Form.Label>Stability</Form.Label>
                <Form.Control
                  {...register('stability')}
                  type="number"
                  onWheel={(e) => e.target.blur()}
                  step="any"
                  min={0}
                  max={1}
                  defaultValue={0.8}
                  disabled={submitted || id}
                ></Form.Control>
              </Form.Group>
              <Form.Group>
                <Form.Label>Minimum Stability</Form.Label>
                <Form.Control
                  {...register('min_stability')}
                  type="number"
                  onWheel={(e) => e.target.blur()}
                  step="any"
                  min={0}
                  max={1}
                  defaultValue={0.2}
                  disabled={submitted || id}
                ></Form.Control>
              </Form.Group>
              <Form.Group>
                <Form.Label>Combined Stability</Form.Label>
                <Form.Control
                  {...register('combined_stability')}
                  type="number"
                  onWheel={(e) => e.target.blur()}
                  step="any"
                  min={0}
                  max={1}
                  defaultValue={1}
                  disabled={submitted || id}
                ></Form.Control>
              </Form.Group>
              <Form.Group>
                <Controller
                  name="allow_stability_drop"
                  control={control}
                  render={({ field }) => (
                    <Form.Check
                      {...field}
                      id="allow_stability_drop"
                      type="checkbox"
                      label={'Allow Stability Drop'}
                      defaultChecked={false}
                      checked={field.checked}
                      disabled={submitted || id}
                    />
                  )}
                />
              </Form.Group>
            </fieldset>

            <fieldset className="border rounded p-2 mb-3">
              <legend className="font-weight-bold">Decomposition</legend>
              <Form.Group>
                <Controller
                  name="make_decomposition_plots"
                  control={control}
                  render={({ field }) => (
                    <Form.Check
                      {...field}
                      id="make_decomposition_plots"
                      type="checkbox"
                      label={'Make Decomposition Plots'}
                      defaultChecked={true}
                      checked={field.checked}
                      disabled={submitted || id}
                    />
                  )}
                />
              </Form.Group>
              <Form.Group>
                <Controller
                  name="collapse_to_SBS96"
                  control={control}
                  render={({ field }) => (
                    <Form.Check
                      {...field}
                      id="collpase_to_SBS96"
                      type="checkbox"
                      label={'Collapse to SBS96'}
                      defaultChecked={true}
                      checked={field.checked}
                      disabled={submitted || id}
                    />
                  )}
                />
              </Form.Group>
            </fieldset>

            <fieldset className="border rounded p-2 mb-3">
              <legend className="font-weight-bold">Others</legend>
              <Form.Group>
                <Controller
                  name="get_all_signature_matrices"
                  control={control}
                  render={({ field }) => (
                    <Form.Check
                      {...field}
                      id="get_all_signature_matrices"
                      type="checkbox"
                      label={'Get All Signature Matrices'}
                      defaultChecked={true}
                      checked={field.checked}
                      disabled={submitted || id}
                    />
                  )}
                />
              </Form.Group>
              <Form.Group>
                <Controller
                  name="export_probabilities"
                  control={control}
                  render={({ field }) => (
                    <Form.Check
                      {...field}
                      id="export_probabilities"
                      type="checkbox"
                      label={'Export Probabilities'}
                      defaultChecked={true}
                      checked={field.checked}
                      disabled={submitted || id}
                    />
                  )}
                />
              </Form.Group>
            </fieldset>
          </div>

          <hr className="mb-3" />
          <div>
            <Form.Group controlId="email">
              <Form.Label>
                Email <span style={{ color: 'crimson' }}>*</span>
              </Form.Label>
              <Controller
                name="email"
                control={control}
                rules={{
                  required: true,
                  pattern:
                    /^[a-zA-Z0-9.!#$%&'*+/=?^_`{|}~-]+@[a-zA-Z0-9-]+(?:\.[a-zA-Z0-9-]+)*$/,
                }}
                render={({ field }) => (
                  <Form.Control
                    {...field}
                    aria-label="email"
                    placeholder="Enter Email"
                    type="email"
                    disabled={submitted || id}
                    isInvalid={errors.email}
                  />
                )}
              />
              <Form.Control.Feedback type="invalid">
                {errors.email && 'Email required'}
              </Form.Control.Feedback>
              <Form.Text className="text-muted">
                <i>You will receive an email when your job is complete</i>
              </Form.Text>
            </Form.Group>
            <Form.Group controlId="jobName">
              <Form.Label>
                Job Name <span style={{ color: 'crimson' }}>*</span>
              </Form.Label>
              <Form.Control
                {...register('jobName', { required: true })}
                placeholder="Enter Job Name"
                disabled={submitted || id}
              />
            </Form.Group>
          </div>
        </div>
        <Row className="mt-3">
          <Col>
            <Button
              disabled={fetchingSeqmatrixOptions}
              className="w-100"
              variant="secondary"
              onClick={() => handleReset()}
            >
              Reset
            </Button>
          </Col>
          <Col>
            <Button
              disabled={
                fetchingSeqmatrixOptions || submitted || id || warning !== ''
              }
              className="w-100"
              variant="primary"
              type="submit"
            >
              Submit
            </Button>
          </Col>
        </Row>
      </Form>
    </div>
  );
}
