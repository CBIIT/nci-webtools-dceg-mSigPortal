import { useState, useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import Select from '../../../controls/select/selectHookForm';
import { useForm, Controller, useFieldArray } from 'react-hook-form';
import { NavHashLink } from 'react-router-hash-link';
import { useProfileComparisonReferenceQuery } from './apiSlice';
import { useSignatureOptionsQuery } from '../../../../services/store/rootApi';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { defaultMatrix } from '../../../../services/utils';
import Plotly from '../../../controls/plotly/plot/plot';
import { useSeqmatrixOptionsQuery } from '../../../../services/store/rootApi';
import { defaultProfile2 } from '../../../../services/utils';

export default function PcReference({ state }) {
  const [params, setParams] = useState(null);
  const { study, cancer, strategy, id, source } = state;

  // query form options
  const { data: seqmatrixOptions } = useSeqmatrixOptionsQuery(
    {
      ...(source == 'public'
        ? { study: study.value, cancer: cancer.value, strategy: strategy.value }
        : { userId: id }),
    },
    { skip: source == 'user' ? !id : !study }
  );
  const { data: signatureOptions, isFetching: fetchingSignatureOptions } =
    useSignatureOptionsQuery();
  // query plot
  const {
    data: plot,
    error,
    isFetching,
  } = useProfileComparisonReferenceQuery(params, {
    skip: !params,
  });

  // define forms
  const defaultValues = {
    profile: '',
    sample: '',
    signatureSet: '',
    compare: [{ proportion: 1, signature: '' }],
  };
  const {
    control,
    handleSubmit,
    watch,
    setValue,
    formState: { errors: formErrors },
  } = useForm({ defaultValues });
  // dynamic signature form
  const { fields, append, remove } = useFieldArray({
    control,
    name: 'compare',
  });
  const { profile, sample, signatureSet, compare } = watch();

  // create form Options
  const profileOptions = seqmatrixOptions
    ? [
        ...new Set(
          seqmatrixOptions
            .map((e) => e.profile)
            .sort((a, b) => b.localeCompare(a))
        ),
      ].map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  const sampleOptions = (profile) =>
    seqmatrixOptions && profile
      ? [
          ...new Set(
            seqmatrixOptions
              .filter((e) => e.profile == profile.value)
              .map((e) => e.sample)
              .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
          ),
        ].map((e) => ({
          label: e,
          value: e,
        }))
      : [];

  const signatureSetOptions = (profile) =>
    signatureOptions && profile
      ? [
          ...new Set(
            signatureOptions
              .filter(
                (e) =>
                  e.profile == profile.value &&
                  e.matrix == defaultMatrix(profile.value, ['96', '78', '83'])
              )
              .map((e) => e.signatureSetName)
          ),
        ]
          .sort((a, b) =>
            a.localeCompare(b, undefined, { sensitivity: 'base' })
          )
          .map((e) => ({
            label: e,
            value: e,
          }))
      : [];

  const signatureNameOptions = (profile, signatureSet) =>
    signatureOptions && profile && signatureSet
      ? [
          ...new Set(
            signatureOptions
              .filter(
                (e) =>
                  e.profile == profile.value &&
                  e.matrix ==
                    defaultMatrix(profile.value, ['96', '78', '83']) &&
                  e.signatureSetName == signatureSet.value
              )
              .map((e) => e.signatureName)
          ),
        ]
          .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
          .map((e) => ({ label: e, value: e }))
      : [];

  // set initial form options
  useEffect(() => {
    if (!profile && profileOptions.length && signatureOptions)
      handleProfile(defaultProfile2(profileOptions));
  }, [profileOptions, signatureOptions]);

  function onSubmit(data) {
    const { profile, sample, signatureSet, compare } = data;
    const signatureName = compare.map((e) => e.signature.value).join(';');
    const scalarValue = compare.map((e) => e.proportion).join(';');

    const spectrumQueryParams = {
      ...(source == 'public' && {
        study: study.value,
        cancer: cancer.value,
        strategy: strategy.value,
      }),
      ...(source == 'user' && { userId: id }),
      profile: profile.value,
      sample: sample.value,
      matrix:
        data.profile.value === 'SBS'
          ? '96'
          : data.profile.value === 'DBS'
          ? '78'
          : '83',
    };

    const signatureQueryParams = {
      profile: profile.value,
      matrix:
        data.profile.value === 'SBS'
          ? '96'
          : data.profile.value === 'DBS'
          ? '78'
          : '83',
      signatureSetName: signatureSet.value,
      signatureName,
      scalarValue,
    };

    setParams({ spectrumQueryParams, signatureQueryParams });
  }

  function handleProfile(profile) {
    const samples = sampleOptions(profile);
    const signatureSets = signatureSetOptions(profile);

    setValue('profile', profile);
    setValue('sample', samples.length ? samples[0] : null);
    handleSignatureSet(profile, signatureSets.length ? signatureSets[0] : null);
  }

  function handleSignatureSet(profile, signatureSet) {
    const signatureOptions = signatureNameOptions(profile, signatureSet);

    setValue('signatureSet', signatureSet);
    fields.forEach((e, i) => {
      setValue(
        `compare.${i}.signature`,
        signatureOptions[i] || signatureOptions[0]
      );
    });
  }

  return (
    <div>
      <p className="p-3 m-0">
        Input a [Profile Type], [Sample Name], [Reference Signature Set], and
        [Compare Signatures] parameter to generate the mutational profile of the
        sample, the signature profile from the reference set, and the difference
        between two mutational profiles. The [Compare Signatures] can be a
        single reference signature name or a combined multiple reference
        signature name with different contributions (see example below [Compare
        Signatures]).
      </p>

      <hr />
      <LoadingOverlay active={isFetching} />
      <Form className="p-3" onSubmit={handleSubmit(onSubmit)}>
        <Form.Row>
          <Col lg="auto">
            <Select
              disabled={!profileOptions.length}
              name="profile"
              label="Profile Type"
              options={profileOptions}
              onChange={handleProfile}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              disabled={!profile}
              name="sample"
              label="Sample Name"
              options={sampleOptions(profile)}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              disabled={fetchingSignatureOptions}
              name="signatureSet"
              label="Reference Signature Set"
              options={signatureSetOptions(profile)}
              onChange={(e) => handleSignatureSet(profile, e)}
              control={control}
            />
          </Col>

          <Col lg="auto" className="d-flex justify-content-end">
            <Button
              className="mt-4 mb-4"
              disabled={
                !profile ||
                !sample ||
                !signatureSet ||
                !compare ||
                fetchingSignatureOptions ||
                fetchingSignatureOptions
              }
              variant="primary"
              type="submit"
            >
              Calculate
            </Button>
          </Col>
        </Form.Row>
        {fields.map((item, index) => (
          <Form.Row key={item.id}>
            <Col sm="1">
              <Form.Group controlId={`proportion-${index}`}>
                <Form.Label>Proportion</Form.Label>
                <Controller
                  name={`compare.${index}.proportion`}
                  rules={{ required: true, min: 0, max: 1 }}
                  control={control}
                  render={({ field }) => (
                    <Form.Control
                      {...field}
                      min={0}
                      max={1}
                      step={0.1}
                      type="number"
                      placeholder="Decimal value between 0 to 1 (e.g. 0.5)"
                    />
                  )}
                />
              </Form.Group>
            </Col>
            <Col sm="auto">
              <Select
                disabled={!signatureNameOptions(profile, signatureSet).length}
                name="signature"
                label="Signature Name"
                value={compare[index].signature}
                options={signatureNameOptions(profile, signatureSet)}
                onChange={(e) => setValue(`compare.${index}.signature`, e)}
                control={control}
              />
            </Col>
            {index > 0 && (
              <Col sm="auto" className="my-auto mx-3">
                <Button
                  className="text-danger p-0"
                  variant="link"
                  title="Remove Signature"
                  onClick={() => remove(index)}
                >
                  - Remove Signature
                </Button>
              </Col>
            )}
            {index == fields.length - 1 && (
              <Col sm="auto" className="my-auto">
                <Button
                  className="p-0"
                  variant="link"
                  title="Add Signature"
                  aria-label="Add Signature"
                  onClick={() =>
                    append({
                      ...defaultValues.compare[0],
                      signature:
                        signatureNameOptions(profile, signatureSet)[
                          index + 1
                        ] || signatureNameOptions(profile, signatureSet)[0],
                    })
                  }
                >
                  + Add Signature
                </Button>
              </Col>
            )}
          </Form.Row>
        ))}
        {sampleOptions(profile).length < 2 && (
          <Row>
            <Col>Unavailable - More than one Sample Required</Col>
          </Row>
        )}
      </Form>
      <div id="pcReferencePlot">
        {error && (
          <>
            <hr />
            <div className="p-3">
              <p>An error has occured. Please verify your input.</p>
              <p>{error.data}</p>
            </div>
          </>
        )}
        {plot && (
          <>
            <hr />
            <Plotly
              className="w-100"
              data={plot.traces}
              layout={plot.layout}
              config={plot.config}
            />
            <div className="p-3">
              <p>
                The plot above shows the mutational profiles of a selected
                sample, the signature from the selected reference signature set,
                and the difference between them. The text at the top of the plot
                indicates the profile similarity calculated using Residual Sum
                of Squares (RSS) and cosine similarity methods.
              </p>
              <p>
                RSS measures the discrepancy between two mutational profiles.
                Cosine similarity measures how similar two mutational profiles
                are. For example, two identical mutational profiles will have
                RSS = 0 and Cosine similarity = 1. For additional information
                about RSS and cosine similarity, click{' '}
                <NavHashLink to="/faq#cosine-similarity">here</NavHashLink>.
              </p>
            </div>
          </>
        )}
      </div>
    </div>
  );
}
