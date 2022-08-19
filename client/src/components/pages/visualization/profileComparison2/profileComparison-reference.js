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
import {
  useProfileComparisonReference1Query,
  useProfileComparisonReference2Query,
  usePcSignatureSetsQuery,
  usePcSignatureNamesQuery,
} from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import SvgContainer from '../../../controls/svgContainer/svgContainer';
import { defaultMatrix } from '../../../../services/utils';

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
  //   calculate plot
  const { data1, error1, isFetching1 } = useProfileComparisonReference1Query(
    calculationQuery,
    { skip: !calculationQuery }
  );
  const { data2, error2, isFetching2 } = useProfileComparisonReference2Query(
    calculationQuery,
    { skip: !calculationQuery }
  );

  const data = { ...data1, ...data2 };
  console.log(data);

  // declare form Options
  const profileOptions = matrixData.length
    ? [...new Set(matrixData.map((e) => e.profile))].map((e) => ({
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
    const params =
      source == 'user'
        ? {
            fn: 'profileComparisonRefSig',
            args: {
              profileType: profile.value,
              sampleName: sample.value,
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
            fn: 'profileComparisonRefSigPublic',
            args: {
              profileType: profile.value,
              sampleName: sample.value,
              signatureSet: signatureSet.value,
              compare: compare,
              study: study.value,
              cancerType: cancer.value,
              experimentalStrategy: strategy.value,
            },
            projectID,
          };
    setCalculationQuery(params);
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
        {data && (
          <>
            <hr />
            {/* <SvgContainer
              className="p-3"
              downloadName={data.output.plotPath.split('/').slice(-1)[0]}
              plotPath={'web/results/' + data.output.plotPath}
              txtPath={`web/results/${data.output.plotPath}`}
            /> */}
            <div className="p-3">
              <p>
                The plot above shows the mutational profiles of two selected
                samples, as well as the difference between them. The text at the
                top of the plot indicates the profile similarity calculated
                using Residual Sum of Squares (RSS) and cosine similarity
                methods.
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
