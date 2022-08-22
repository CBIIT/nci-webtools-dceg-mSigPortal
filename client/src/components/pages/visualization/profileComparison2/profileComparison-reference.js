import { useState, useEffect } from 'react';
import {
  Form,
  Row,
  Col,
  Button,
  Popover,
  OverlayTrigger,
} from 'react-bootstrap';
import Plot from 'react-plotly.js';
import Select from '../../../controls/select/selectForm';
import { useForm, Controller } from 'react-hook-form';
import { NavHashLink } from 'react-router-hash-link';
import { useSelector, useDispatch } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faInfoCircle } from '@fortawesome/free-solid-svg-icons';
import { actions } from '../../../../services/store/visualization';
import {
  useProfileComparisonReference1Query,
  useProfileComparisonReference2Query,
  usePcSignatureSetsQuery,
  usePcSignatureNamesQuery,
} from './apiSlice';
import { cloneDeep } from 'lodash';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import SvgContainer from '../../../controls/svgContainer/svgContainer';
import { defaultMatrix } from '../../../../services/utils';
import pcReferenceSignature_SBS from '../../../controls/plotly/profileComparision/pcReferenceSignature_SBS';
import pcBetweenSamples_DBS from '../../../controls/plotly/profileComparision/pcBetweenSamples_DBS';
import pcBetweenSamples_ID from '../../../controls/plotly/profileComparision/pcBetweenSamples_ID';

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
  const { source, matrixData, svgList, matrixList, projectID } = store.main;
  const { referenceForm } = store.profileComparison;
  console.log('study:');
  console.log(study);
  console.log('cancer:');
  console.log(cancer);
  console.log('strategy');
  console.log(strategy);
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

  const [calculationQuery1, setCalculationQuery1] = useState('');
  const [calculationQuery2, setCalculationQuery2] = useState('');
  const [signatureSetQuery, setSignatureSetQuery] = useState('');
  const [signatureNamesQuery, setSignatureNamesQuery] = useState('');

  // get signature sets
  const { data: signatureSetOptions, isFetching: fetchingSigSets } =
    usePcSignatureSetsQuery(signatureSetQuery, {
      skip: !signatureSetQuery,
    });
  // get signature names in set
  const {
    data: signatureNameOptions,
    isFetching: fetchingSigNames,
    refetch: refetchSignatureNames,
  } = usePcSignatureNamesQuery(signatureNamesQuery, {
    skip: !signatureNamesQuery,
  });
  //   seqmatrix api
  const {
    data: data1,
    error: error1,
    isFetching: isFetching1,
  } = useProfileComparisonReference1Query(calculationQuery1, {
    skip: !calculationQuery1,
  });

  //signature api
  const {
    data: data2,
    error: error2,
    isFetching: isFetching2,
  } = useProfileComparisonReference2Query(calculationQuery2, {
    skip: !calculationQuery2,
  });
  let datatotal;
  if ((data1 != null) & (data2 != null)) {
    console.log('data1');
    console.log(data1);
    console.log('data2');
    console.log(data2);
    console.log(data1.arg);
    console.log(data2.arg);
    const compares = [data1.arg.sample, data2.arg.signatureName];
    console.log(compares);
    const dat1temp = Object.values(data1);
    const sample1 = dat1temp.slice(0, dat1temp.length - 1);
    console.log(sample1);
    const dat2temp = Object.values(data2);
    const sample2 = dat2temp.slice(0, dat2temp.length - 1);
    console.log(sample2);
    datatotal = [sample1, sample2];
    if (data1.arg.profile === 'SBS') {
      pcReferenceSignature_SBS(compares, sample1, sample2, 'signatures');
    } else if (data1.arg.profile === 'DBS') {
      pcBetweenSamples_DBS(compares, sample1, sample2, 'signatures');
    } else {
      pcBetweenSamples_ID(compares, sample1, sample2, 'signatures');
    }
  }

  // declare form Options
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

  const sampleOptions = getSampleOptions(profile);

  // set inital profile
  useEffect(() => {
    if (!profile && profileOptions.length) handleProfile(profileOptions[0]);
  }, [profileOptions]);
  // set intital signature set
  useEffect(() => {
    if (signatureSetOptions) setValue('signatureSet', signatureSetOptions[0]);
  }, [signatureSetOptions]);
  // set initial signature
  useEffect(() => {
    if (signatureNameOptions) setValue('compare', signatureNameOptions[0]);
  }, [signatureNameOptions]);
  // get signature sets when profile is selected
  useEffect(() => {
    if (profile) {
      setSignatureSetQuery({
        profile: profile.value,
        matrix: defaultMatrix(profile.value, ['96', '78', '83']),
      });
    }
  }, [profile]);
  // get signature names when signature set is selected
  useEffect(() => {
    if (signatureSet) {
      setSignatureNamesQuery({
        profile: profile.value,
        matrix: defaultMatrix(profile.value, ['96', '78', '83']),
        signatureSetName: signatureSet.value,
      });
    }
  }, [signatureSet]);

  // get samples filtered by selected profile
  function getSampleOptions(profile) {
    return matrixData && profile
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
  }

  function onSubmit(data) {
    mergeForm(data);

    const { profile, sample, signatureSet, compare } = data;
    const params1 =
      source == 'user'
        ? {
            fn: 'profileComparisonRefSig',
            args: {
              study: study.value,
              cancer: cancer.value,
              strategy: strategy.value,
              profile: profile.value,
              sample: sample.value,
              signatureSet: signatureSet.value,
              compare: compare,
              matrixFile: matrixList.filter(
                (e) =>
                  e.profile == profile.value &&
                  e.matrix == defaultMatrix(profile.value, ['96', '78', '83'])
              )[0].Path,
            },
            projectID,
          }
        : {
            study: study.value,
            cancer: cancer.value,
            strategy: strategy.value,
            profile: profile.value,
            sample: sample.value,
            matrix:
              data.profile.value === 'SBS'
                ? '96'
                : data.profile.value === 'DBS'
                ? '78'
                : '83',
          };
    const params2 =
      source == 'user'
        ? {
            fn: 'profileComparisonRefSig',
            args: {
              strategy: strategy.value,
              profileType: profile.value,
              sample: sample.value,
              signatureSetName: signatureSet.value,
              signatureName: compare,
              matrixFile: matrixList.filter(
                (e) =>
                  e.profile == profile.value &&
                  e.matrix == defaultMatrix(profile.value, ['96', '78', '83'])
              )[0].Path,
            },
            projectID,
          }
        : {
            strategy: strategy.value,
            profile: profile.value,
            matrix:
              data.profile.value === 'SBS'
                ? '96'
                : data.profile.value === 'DBS'
                ? '78'
                : '83',
            signatureSetName: signatureSet.value,
            signatureName: compare,
          };
    setCalculationQuery1(params1);
    setCalculationQuery2(params2);
  }

  function handleProfile(e) {
    const samples = getSampleOptions(e);

    setValue('profile', e);
    if (samples.length) {
      setValue('sample', samples[0]);
    }
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
          {signatureNameOptions ? (
            signatureNameOptions
              .filter((e) =>
                search ? e.toLowerCase().includes(search.toLowerCase()) : true
              )
              .map((signature) => (
                <Button
                  key={signature}
                  variant="link"
                  onClick={() => {
                    if (compare.length) {
                      setValue('compare', compare + `;1*${signature}`);
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
      <LoadingOverlay active={isFetching1} />
      <Form className="p-3" onSubmit={handleSubmit(onSubmit)}>
        <Row>
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
              disabled={sampleOptions.length < 2}
              name="sample"
              label="Sample Name"
              options={sampleOptions}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              disabled={sampleOptions.length < 2 || fetchingSigSets}
              name="signatureSet"
              label="Reference Signature Set"
              options={signatureSetOptions}
              onChange={(e) => {
                setValue('signatureSet', e);
                resetField('compare');
                refetchSignatureNames();
              }}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Form.Group controlId="compare">
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
                    disabled={fetchingSigNames}
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
          <Col lg="auto" className="d-flex">
            <Button
              className="mt-auto mb-3"
              disabled={!profile || !sample || !signatureSet || !compare}
              variant="primary"
              type="submit"
            >
              Calculate
            </Button>
          </Col>
        </Row>
        {sampleOptions.length < 2 && (
          <Row>
            <Col>Unavailable - More than one Sample Required</Col>
          </Row>
        )}
      </Form>
      <div id="pcReferencePlot">
        {error1 && (
          <>
            <hr />
            <div className="p-3">
              <p>An error has occured. Please verify your input.</p>
              <p>{error1.data}</p>
            </div>
          </>
        )}
        {data1 && data2 && (
          <>
            <hr />

            <Plot
              className="w-100"
              data={cloneDeep([data1, data2].traces)}
              layout={cloneDeep([data1, data2].layout)}
              config={cloneDeep([data1, data2].config)}
              useResizeHandler
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
