import { useState, useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import Select from '../../../controls/select/selectForm';
import { useForm } from 'react-hook-form';
import { NavHashLink } from 'react-router-hash-link';
import { useSelector, useDispatch } from 'react-redux';
import { actions } from '../../../../services/store/visualization';
import { useProfileComparisonPublicQuery } from './apiSlice';
import { useSeqmatrixOptionsQuery } from '../../../../services/store/rootApi';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plotly from '../../../controls/plotly/plot/plot';

export default function PcPublic() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.visualization);
  const mergeForm = (state) =>
    dispatch(
      actions.mergeVisualization({
        profileComparison: { ...store.profileComparison, publicForm: state },
      })
    );

  const { matrixData, id } = store.main;
  const { publicForm } = store.profileComparison;

  const [params, setParams] = useState(null);

  const { data, error, isFetching } = useProfileComparisonPublicQuery(params, {
    skip: !params,
  });

  const { data: publicOptions, isFetching: fetchingPublicOptions } =
    useSeqmatrixOptionsQuery();

  const studyOptions = publicOptions
    ? [...new Set(publicOptions.map((e) => e.study))].sort().map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  const cancerOptions = (study) => {
    if (publicOptions && study?.value) {
      return [
        ...[
          ...new Set(
            publicOptions
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
    } else {
      return [];
    }
  };

  function publicSampleOptions(study, cancer) {
    return publicOptions && study && cancer
      ? [
          ...new Set(
            publicOptions
              .filter(
                (e) => e.study === study.value && e.cancer === cancer.value
              )
              .map((e) => e.sample)
              .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
          ),
        ].map((e) => ({
          label: e,
          value: e,
        }))
      : [];
  }

  const { control, handleSubmit, watch, setValue } = useForm({
    defaultValues: publicForm,
  });

  const { profile, matrix, userSample, study, cancer, publicSample } = watch();

  const supportedProfiles = ['SBS', 'DBS', 'ID'];
  const supportedMatrices = [96, 192, 78, 83];

  const profileOptions = matrixData.length
    ? [
        ...new Set(
          matrixData.map((e) => e.profile).sort((a, b) => b.localeCompare(a))
        ),
      ]
        .filter((e) => supportedProfiles.includes(e))
        .map((e) => ({
          label: e,
          value: e,
        }))
    : [];

  const matrixOptions = (profile) =>
    matrixData.length && profile
      ? [
          ...new Set(
            matrixData
              .filter((e) => e.profile == profile.value)
              .map((e) => e.matrix)
          ),
        ]
          .filter((e) => supportedMatrices.includes(e))
          .sort((a, b) => a - b)
          .map((e) => ({
            label: e,
            value: e,
          }))
      : [];

  const userSampleOptions = (profile, matrix) =>
    matrixData.length && profile
      ? [
          ...new Set(
            matrixData
              .filter(
                (e) => e.profile == profile.value && e.matrix == matrix.value
              )
              .map((e) => e.sample)
          ),
        ].map((e) => ({
          label: e,
          value: e,
        }))
      : [];

  // set inital parameters
  useEffect(() => {
    if (!profile && profileOptions.length) handleProfile(profileOptions[0]);
  }, [profileOptions]);
  useEffect(() => {
    if (!study && studyOptions.length) handleStudy(studyOptions[0]);
  }, [studyOptions]);

  function handleProfile(profile) {
    const matrices = matrixOptions(profile);
    setValue('profile', profile);
    handleMatrix(profile, matrices[0]);
  }

  function handleMatrix(profile, matrix) {
    const samples = userSampleOptions(profile, matrix);
    setValue('matrix', matrix);
    setValue('userSample', samples[0]);
  }

  function handleStudy(study) {
    const cancers = cancerOptions(study);
    setValue('study', study);
    handleCancer(study, cancers[0]);
  }

  function handleCancer(study, cancer) {
    const samples = publicSampleOptions(study, cancer);
    setValue('cancer', cancer);
    setValue('publicSample', samples[0]);
  }

  function onSubmit(data) {
    const userParams = {
      profile: data.profile.value,
      matrix: data.matrix.value,
      sample: data.userSample.value,
      userId: id,
    };
    const publicParams = {
      profile: data.profile.value,
      matrix: data.matrix.value,
      study: data.study.value,
      cancer: data.cancer.value,
      sample: data.publicSample.value,
    };
    setParams({ userParams, publicParams });
    mergeForm(data);
  }

  return (
    <div>
      <p className="p-3 m-0">
        Input a [Profile Type], [Matrix Size], [Sample Name] (from your input
        data), [Study], [Cancer Type], and [Public Sample Name] to generate the
        mutational profile of the input sample, the sample from the selected
        public data, and the difference between them.
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
              options={profileOptions}
              onChange={handleProfile}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              disabled={isFetching}
              name="matrix"
              label="Matrix size"
              options={matrixOptions(profile)}
              onChange={(e) => handleMatrix(profile, e)}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              disabled={!userSampleOptions(profile, matrix)}
              name="userSample"
              label="Sample"
              options={userSampleOptions(profile, matrix)}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="study"
              label="Study"
              options={studyOptions}
              onChange={handleStudy}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="cancer"
              label="Cancer Type or Gorup"
              options={cancerOptions(study)}
              onChange={(e) => handleCancer(study, e)}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              disabled={!publicSampleOptions(study, cancer)}
              name="publicSample"
              label="Sample Name"
              options={publicSampleOptions(study, cancer)}
              control={control}
            />
          </Col>
          <Col lg="auto" className="d-flex">
            <Button
              className="mt-auto mb-3"
              disabled={
                !profile ||
                !matrix ||
                !userSample ||
                !study ||
                !cancer ||
                !publicSample
              }
              variant="primary"
              type="submit"
            >
              Calculate
            </Button>
          </Col>
        </Row>
      </Form>
      <div id="plot">
        {error && (
          <>
            <hr />
            <div className="p-3">
              <p>An error has occured. Please verify your input.</p>
              <p>{error.data || error.message || error}</p>
            </div>
          </>
        )}
        {data && (
          <>
            <hr />
            <Plotly
              className="w-100"
              data={data.traces}
              layout={data.layout}
              config={data.config}
            />
            <div className="p-3">
              <p>
                The plot above shows the mutational profile of the input sample,
                the sample from the selected public data, and the difference
                between them. The text on the top of the plot indicates the
                profile similarity calculated using Residual Sum of Squares
                (RSS) and cosine similarity methods.
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
