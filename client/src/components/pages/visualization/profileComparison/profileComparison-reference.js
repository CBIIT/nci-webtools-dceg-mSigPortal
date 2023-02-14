import { useState, useEffect } from 'react';
import {
  Form,
  Row,
  Col,
  Button,
  Popover,
  OverlayTrigger,
} from 'react-bootstrap';
import Select from '../../../controls/select/selectForm';
import { useForm, Controller } from 'react-hook-form';
import { NavHashLink } from 'react-router-hash-link';
import { useSelector, useDispatch } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faInfoCircle } from '@fortawesome/free-solid-svg-icons';
import { actions } from '../../../../services/store/visualization';
import { useProfileComparisonReferenceQuery } from './apiSlice';
import { useSignatureOptionsQuery } from '../../../../services/store/rootApi';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { defaultMatrix } from '../../../../services/utils';
import Plotly from '../../../controls/plotly/plot/plot';

export default function PcReference() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.visualization);
  const mergeForm = (state) =>
    dispatch(
      actions.mergeVisualization({
        profileComparison: { ...store.profileComparison, referenceForm: state },
      })
    );

  const { study, cancer, strategy } = store.publicForm;
  const { source, matrixData, id } = store.main;
  const { referenceForm } = store.profileComparison;

  // main form
  const {
    control,
    handleSubmit,
    watch,
    setValue,
    resetField,
    formState: { errors: formErrors },
  } = useForm({
    defaultValues: referenceForm,
  });
  const { profile, sample, signatureSet, compare } = watch();

  //   signature search form
  const { control: searchControl, watch: watchSearch } = useForm({
    defaultValues: { search: '' },
  });
  const { search } = watchSearch();

  const [calculationQuery, setCalculationQuery] = useState('');

  // get signature options
  const { data: signatureOptions, isFetching: fetchingSignatureOptions } =
    useSignatureOptionsQuery();

  // get plot data
  const {
    data: plot,
    error,
    isFetching,
  } = useProfileComparisonReferenceQuery(calculationQuery, {
    skip: !calculationQuery,
  });

  // create form Options
  const profileOptions = matrixData.length
    ? [
        ...new Set(
          matrixData.map((e) => e.profile).sort((a, b) => b.localeCompare(a))
        ),
      ].map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  const sampleOptions = (profile) =>
    matrixData && profile
      ? [
          ...new Set(
            matrixData
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
        ].sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
      : [];

  // set inital form options
  useEffect(() => {
    if (!profile && profileOptions.length && signatureOptions)
      handleProfile(profileOptions[0]);
  }, [profileOptions, signatureOptions]);

  function onSubmit(data) {
    mergeForm(data);
    const { profile, sample, signatureSet, compare } = data;
    const params_spectrum = {
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
    const paramsArray = compare.split(';');
    let paramsScalar = [];
    let paramsSig = [];
    for (var i = 0; i < paramsArray.length; i++) {
      if (paramsArray[i].includes('*')) {
        const parts = paramsArray[i].split('*');
        paramsScalar.push(parts[0]);
        paramsSig.push(parts[1]);
      } else {
        paramsScalar.push('1');
        paramsSig.push(paramsArray[i]);
      }
    }

    const params_signature = {
      profile: profile.value,
      matrix:
        data.profile.value === 'SBS'
          ? '96'
          : data.profile.value === 'DBS'
          ? '78'
          : '83',
      signatureSetName: signatureSet.value,
      signatureName: paramsSig.join(';'),
      scalarValue: paramsScalar.join(';'),
    };
    const params_signature_scalar = {
      arrayScalar: paramsScalar,
      arraySignature: paramsSig,
      paramsSignatureScalar: compare,
    };

    setCalculationQuery({
      params_spectrum,
      params_signature,
      params_signature_scalar,
    });
  }

  function handleProfile(profile) {
    const samples = sampleOptions(profile);
    const signatureSets = signatureSetOptions(profile);

    setValue('profile', profile);
    setValue('sample', samples.length ? samples[0] : null);
    handleSignatureSet(profile, signatureSets.length ? signatureSets[0] : null);
  }

  function handleSignatureSet(profile, signatureSet) {
    const names = signatureNameOptions(profile, signatureSet);

    setValue('signatureSet', signatureSet);
    setValue('compare', names[0]);
  }

  const signatureListPopover = signatureSet && (
    <Popover id="popover-basic" style={{ minWidth: '400px' }}>
      <Popover.Title as="h3">{signatureSet.value}</Popover.Title>
      <Popover.Content>
        <Controller
          name="search"
          control={searchControl}
          render={({ field }) => (
            <Form.Control {...field} placeholder="Search Signatures" />
          )}
        />
        <div>
          {signatureNameOptions(profile, signatureSet) ? (
            signatureNameOptions(profile, signatureSet)
              .filter((e) =>
                search ? e.toLowerCase().includes(search.toLowerCase()) : true
              )
              .map((signature) => (
                <Button
                  key={signature}
                  variant="link"
                  onClick={() => {
                    if (compare.length) {
                      setValue('compare', compare + `;${signature}`);
                    } else {
                      setValue('compare', signature);
                    }
                  }}
                >
                  {signature}
                </Button>
              ))
          ) : (
            <p>No Signatures Available</p>
          )}
        </div>
      </Popover.Content>
    </Popover>
  );

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
        <Row>
          <Col lg="auto">
            <Select
              disabled={!profileOptions.length}
              name="profile"
              label="Profile Type"
              className="m-auto"
              options={profileOptions}
              onChange={handleProfile}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              disabled={!profile}
              className="m-auto"
              name="sample"
              label="Sample Name"
              options={sampleOptions(profile)}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              disabled={fetchingSignatureOptions}
              className="m-auto"
              name="signatureSet"
              label="Reference Signature Set"
              options={signatureSetOptions(profile)}
              onChange={(e) => handleSignatureSet(profile, e)}
              control={control}
            />
          </Col>
          <Col lg="auto" className="d-flex">
            <Form.Group controlId="compare" className="m-auto">
              <Form.Label>
                Compare Signatures{' '}
                <OverlayTrigger
                  trigger="click"
                  placement="top"
                  overlay={signatureListPopover}
                  rootClose
                >
                  <Button
                    disabled={!signatureSet}
                    aria-label="compare signatures info"
                    variant="link"
                    className="p-0 font-weight-bold "
                  >
                    <FontAwesomeIcon
                      icon={faInfoCircle}
                      style={{ verticalAlign: 'baseline' }}
                    />
                  </Button>
                </OverlayTrigger>
              </Form.Label>
              <Controller
                name="compare"
                control={control}
                rules={{ required: true }}
                render={({ field }) => (
                  <Form.Control
                    {...field}
                    isInvalid={formErrors?.compare}
                    disabled={fetchingSignatureOptions}
                  />
                )}
              />
              <Form.Text className="text-muted">
                (Ex. 0.8*SBS5;0.1*SBS1)
              </Form.Text>
              <Form.Control.Feedback type="invalid">
                Enter a valid signature. Click info icon for options.
              </Form.Control.Feedback>
            </Form.Group>
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
        </Row>
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
