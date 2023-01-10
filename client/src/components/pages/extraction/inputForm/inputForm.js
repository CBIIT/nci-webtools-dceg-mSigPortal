import { useState, useEffect } from 'react';
import { Form, Button, Row, Col } from 'react-bootstrap';
import { useForm, Controller } from 'react-hook-form';
import SelectForm from '../../../controls/select/selectForm';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { useSelector, useDispatch } from 'react-redux';
import { actions as extractionActions } from '../../../../services/store/extraction';
import { actions as modalActions } from '../../../../services/store/modal';
import {
  resetExtractionApi,
  useSeqmatrixOptionsQuery,
  useSignatureOptionsQuery,
  useRefGenomeQuery,
} from '../../../../services/store/rootApi';
import { useUploadMutation, useSubmitMutation } from './apiSlice';

const actions = { ...extractionActions, ...modalActions };

export default function InputForm() {
  const store = useSelector((state) => state.extraction);
  const { submitted } = store.main;
  const { inputFilename, ...inputForm } = store.inputForm;

  const dispatch = useDispatch();
  const mergeForm = (state) =>
    dispatch(actions.mergeExtraction({ inputForm: state }));
  const resetExtraction = (_) => dispatch(actions.resetExtraction());
  const mergeMain = (state) =>
    dispatch(actions.mergeExtraction({ main: state }));

  // query options to populate form
  const {
    data: seqmatrixOptions,
    error: seqmatrixError,
    isFetching: fetchingSeqmatrixOptions,
  } = useSeqmatrixOptionsQuery();
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

  const signatureSetOptions = signatureOptions
    ? [...new Set(signatureOptions.map((e) => e.signatureSetName))].map(
        (e) => ({ label: e, value: e })
      )
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
            ),
          ].map((e) => ({ label: e, value: e })),
        ]
      : [];

  const contextTypeOptions = signatureOptions
    ? [
        { label: 'default', value: 'SBS96' },
        ...[...new Set(signatureOptions.map((e) => e.profile + e.matrix))]
          .sort((a, b) => a.localeCompare(b, 'en', { numeric: true }))
          .map((e) => ({ label: e, value: e })),
      ]
    : [];

  const dataTypeOptions = ['vcf', 'matrix'].map((e) => ({
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
    source: 'user',
    study: { label: 'PCAWG', value: 'PCAWG' },
    cancer: { label: 'Lung-AdenoCA', value: 'Lung-AdenoCA' },
    strategy: { label: 'WGS', value: 'WGS' },
    input_type: { label: 'vcf', value: 'vcf' },
    reference_genome: { label: 'GRCh37', value: 'GRCh37' },
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

    matrix_normalization: 'gmm',
    nmf_init: 'random',
    precision: 'single',
  };

  const sample1 = {
    source: 'user',
    study: { label: 'PCAWG', value: 'PCAWG' },
    cancer: { label: 'Lung-AdenoCA', value: 'Lung-AdenoCA' },
    strategy: { label: 'WGS', value: 'WGS' },
    input_type: { label: 'matrix', value: 'matrix' },
    reference_genome: { label: 'GRCh37', value: 'GRCh37' },
    exome: false,
    context_type: { label: 'SBS96', value: 'SBS96' },
    signatureSetName: {
      label: 'COSMIC_v3.3_Signatures_GRCh37_SBS96',
      value: 'COSMIC_v3.3_Signatures_GRCh37_SBS96',
    },
    signatureName: [{ label: 'all', value: 'all' }],
    extractTool: {
      label: 'SigProfilerExtractor',
      value: 'SigProfilerExtractor',
    },

    minimum_signatures: 1,
    maximum_signatures: 2,
    min_nmf_iterations: 2,
    max_nmf_iterations: 4,
    nmf_test_conv: 2,
    nmf_replicates: 10,

    email: 'test.email@nih.gov',
  };
  const sample2 = {};

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
    exome,
    context_type,
    signatureSetName,
    minSignatures,
    input_type,
  } = watch();

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
    window.location.hash = '#/extraction';
    resetForm(defaultValues);
    resetExtraction();
    dispatch(resetExtractionApi);
  }

  async function onSubmit(data) {
    console.log(data);

    const args = {
      ...(source == 'user' && {
        input_type: data.input_type.value,
        input_data: data.inputFile.name,
      }),
      reference_genome: data.reference_genome.value,
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
    };
    const signatureQuery = {
      signatureSetName: data.signatureSetName.value,
      profile: data.context_type.value.match(/^\D*/)[0],
      matrix: data.context_type.value.match(/\d*$/)[0],
      ...(data.signatureName[0].value != 'all' && {
        signatureName: data.signatureName.map((e) => e.value).join(';'),
      }),
    };

    const formData = new FormData();
    formData.append('inputFile', data.inputFile);

    const { id } = await uploadFiles(formData).unwrap();
    const params = { args, signatureQuery, id, email: data.email };
    const res = await submitForm(params).unwrap();

    console.log(res);

    mergeForm(data);
    mergeMain({ id });
    // try {
    //   mergeMain({ submitted: true, loading: { active: true } });
    //   mergeState(data);
    //   const params = {
    //     study: data.study.value,
    //     cancer: data.cancer.value,
    //     strategy: data.strategy.value,
    //   };

    //   // let matrixData = [];
    //   // for await (const data of paginateQuery(fetchMatrix, params)) {
    //   //   matrixData = [...matrixData, ...data];
    //   // }
    //   const matrixData = []; //await fetchMatrix(params).unwrap();

    //   mergeMain({ matrixData, projectID: crypto.randomUUID() });
    // } catch (error) {
    //   console.log(error);
    //   if (error.originalStatus == 504) {
    //     mergeMain({
    //       error: 'Please Reset Your Parameters and Try again.',
    //     });
    //     mergeError({
    //       visible: true,
    //       message:
    //         'Your submission has timed out. Please try again by submitting this job to a queue instead.',
    //     });
    //   } else {
    //     mergeMain({
    //       error: 'Please Reset Your Parameters and Try again.',
    //     });
    //     mergeError(error.data);
    //   }
    // }
    // mergeMain({ loading: { active: false } });
  }

  return (
    <Form
      onSubmit={handleSubmit(onSubmit)}
      style={{ maxHeight: '900px', overflow: 'hidden auto' }}
    >
      <LoadingOverlay active={fetchingSeqmatrixOptions} />
      {(seqmatrixError || signatureError) && (
        <p>There was an error retrieving public data options</p>
      )}
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
                checked={field.value == 'public'}
                disabled={true}
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
                value={'vcf'}
                checked={field.value == 'user'}
                disabled={submitted}
              />
            )}
          />
        </Form.Group>
        {source == 'public' ? (
          <div>
            <SelectForm
              className="mb-2"
              name="study"
              label="Study"
              disabled={submitted || fetchingSeqmatrixOptions}
              options={studyOptions}
              control={control}
              onChange={handleStudyChange}
            />
            <SelectForm
              className="mb-2"
              name="cancer"
              label="Cancer Type or Group"
              disabled={submitted || fetchingSeqmatrixOptions}
              options={cancerOptions(study)}
              control={control}
              onChange={handleCancerChange}
            />
            <SelectForm
              className="mb-2"
              name="strategy"
              label="Experimental Strategy"
              disabled={submitted || fetchingSeqmatrixOptions}
              options={strategyOptions(study, cancer)}
              control={control}
            />
          </div>
        ) : (
          <div>
            <SelectForm
              name="input_type"
              label="Data Type"
              disabled={submitted}
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
                    disabled={submitted}
                    id="inputFile"
                    title={inputFile?.name || 'Upload Data File...'}
                    label={
                      inputFile?.name || inputFilename || 'Upload Data File...'
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
              onClick={async () => {
                resetForm(sample1);
                const file = 'extraction_sample_SBS96.all';
                const path = 'assets/exampleInput/' + file;
                setValue(
                  'inputFile',
                  new File([await (await fetch(path)).blob()], file)
                );
              }}
            >
              Load Sample
            </Button>
            {/* <Button variant="link" onClick={() => resetForm(sample2)}>
              Sample 2
            </Button> */}
          </div>
        )}
      </div>

      <SelectForm
        name="reference_genome"
        label="Reference Genome Build"
        disabled={submitted}
        options={genomeOptions}
        control={control}
        rules={{ required: input_type == 'vcf' }}
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
              disabled={submitted}
            />
          )}
        />
      </Form.Group>
      <SelectForm
        name="context_type"
        label="Context Type"
        disabled={submitted}
        defaultValue={contextTypeOptions[0]}
        options={contextTypeOptions}
        control={control}
      />
      <SelectForm
        name="signatureSetName"
        label="Reference Signature Set"
        disabled={submitted}
        options={signatureSetOptions}
        control={control}
      />
      <SelectForm
        name="signatureName"
        label="SelectForm Signature Names"
        disabled={submitted}
        options={signatureNameOptions(signatureSetName)}
        control={control}
        onChange={(values, e) => {
          // remove "all" option if a specific signature is selected
          if (e.option.value !== 'all') {
            setValue(
              'signatureName',
              values.filter((e) => e.value !== 'all')
            );
          } else if (e.option.value === 'all') {
            setValue('signatureName', [e.option]);
          } else setValue('signatureName', values);
        }}
        isMulti
      />
      <SelectForm
        name="extractTool"
        label="Extract Tool"
        disabled={submitted}
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
        {/* <fieldset className="border rounded p-2 mb-3">
          <legend className="font-weight-bold">Execution</legend>
          <Form.Group controlId="batch_size">
            <Form.Label>Batch Size</Form.Label>
            <Form.Control
              {...register('batch_size')}
              type="number"
              defaultValue={1}
            />
          </Form.Group>
        </fieldset> */}
        <fieldset className="border rounded p-2 mb-3">
          <legend className="font-weight-bold">NMF Replicates</legend>
          <Form.Group controlId="minSignatures">
            <Form.Label>Minimum Signatures</Form.Label>
            <Form.Control
              {...register('minimum_signatures')}
              type="number"
              min="1"
              max="24"
              defaultValue={1}
            />
          </Form.Group>
          <Form.Group controlId="maxSignatures">
            <Form.Label>Maximum Signatures</Form.Label>
            <Form.Control
              {...register('maximum_signatures')}
              type="number"
              min="1"
              max="25"
              defaultValue={15}
            />
          </Form.Group>
          <Form.Group controlId="nmf_replicates">
            <Form.Label>NMF Replicates Size</Form.Label>
            <Form.Control
              {...register('nmf_replicates')}
              type="number"
              min="1"
              max="1000"
              defaultValue={100}
            />
          </Form.Group>
          <Form.Group controlId="resample">
            <Form.Check
              {...register('resample')}
              label="Resample"
              id="resample"
              defaultChecked={true}
            />
          </Form.Group>
          <Form.Group controlId="seeds">
            <Form.Label>Seeds</Form.Label>
            <Form.Control {...register('seeds')} defaultValue="random" />
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
                <Form.Control {...field} as="select" disabled={submitted}>
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
                <Form.Control {...field} as="select" disabled={submitted}>
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
                <Form.Control {...field} as="select" disabled={submitted}>
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
              {...register('min_nmf_iterations')}
              type="number"
              defaultValue={10000}
            />
          </Form.Group>
          <Form.Group controlId="maxNmfIterations">
            <Form.Label>Maximum NMF Iterations</Form.Label>
            <Form.Control
              {...register('max_nmf_iterations')}
              type="number"
              defaultValue={1000000}
            />
          </Form.Group>
          <Form.Group controlId="nmfTestConv">
            <Form.Label>NMF Test Convergence</Form.Label>
            <Form.Control
              {...register('nmf_test_conv')}
              type="number"
              defaultValue={10000}
            />
          </Form.Group>
          <Form.Group controlId="nmfTolerance">
            <Form.Label>NMF Tolerance</Form.Label>
            <Form.Control
              {...register('nmf_tolerance')}
              type="number"
              defaultValue={1.0e-15}
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
                  disabled={submitted}
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
                  disabled={submitted}
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
                  disabled={submitted}
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
                  disabled={submitted}
                />
              )}
            />
          </Form.Group>
        </fieldset>
      </div>

      <hr className="mb-3" />
      <div>
        <Form.Group controlId="toggleQueue" className="d-flex">
          <Controller
            name="useQueue"
            control={control}
            render={({ field }) => (
              <Form.Check
                {...field}
                id="useQueue"
                type="checkbox"
                label="Long Running Job"
                checked={true}
                disabled={true}
              />
            )}
          />
        </Form.Group>
        <Form.Group controlId="email">
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
                aria-label="Enter Email"
                placeholder="Enter Email"
                size="sm"
                type="email"
                disabled={submitted}
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
      </div>

      <Row>
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
            disabled={fetchingSeqmatrixOptions || submitted}
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
