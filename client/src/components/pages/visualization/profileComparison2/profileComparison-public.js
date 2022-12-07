import { useState, useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import Select from '../../../controls/select/selectForm';
import { useForm } from 'react-hook-form';
import { NavHashLink } from 'react-router-hash-link';
import { useSelector, useDispatch } from 'react-redux';
import { actions } from '../../../../services/store/visualization';
import { useProfileComparisonWithinQuery } from './apiSlice';
import { useSeqmatrixOptionsQuery } from '../../../../services/store/rootApi';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plotly from '../../../controls/plotly/plot/plot';

export default function PcPublic() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.visualization);
  const mergeForm = (state) =>
    dispatch(
      actions.mergeVisualization({
        profileComparison: { ...store.profileComparison, withinForm: state },
      })
    );

  const { matrixData, projectID } = store.main;
  const { withinForm } = store.profileComparison;

  const [params, setParams] = useState(null);

  const { data, error, isFetching } = useProfileComparisonWithinQuery(params, {
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

  function sampleOptions(profile) {
    return matrixData && profile
      ? [
          ...new Set(
            matrixData
              .filter((e) => e.profile === profile.value)
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
    defaultValues: withinForm,
  });

  const { profile, matrix, sample, study, cancer, publicSample } = watch();

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

  // set inital parameters
  useEffect(() => {
    if (!profile && profileOptions.length) handleProfile(profileOptions[0]);
  }, [profileOptions]);

  function onSubmit(data) {
    mergeForm(data);
    const params = {
      userParams: {
        userId: projectID,
        profile: data.profile.value,
        matrix:
          data.profile.value === 'SBS'
            ? '96'
            : data.profile.value === 'DBS'
            ? '78'
            : '83',
        sample: data.sample.value,
      },
      publicParams: {
        study: data.study.value,
        cancer: data.cancer.value,
        sample: data.publicSample.value,
      },
    };
    setParams(params);
  }

  function handleProfile(e) {
    const samples = sampleOptions(e);

    setValue('profile', e);
    if (samples.length) {
      setValue('sample1', samples[0]);
      samples.length > 1
        ? setValue('sample2', samples[1])
        : setValue('sample2', samples[0]);
    }
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
              disabled={true}
              name="matrix"
              label="Sample Matrix size"
              options={[]}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              disabled={!sampleOptions(profile).length}
              name="sample"
              label="Sample Name"
              options={sampleOptions(profile)}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="study"
              label="Study"
              options={studyOptions}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="cancer"
              label="Cancer Type or Gorup"
              options={cancerOptions(study)}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              disabled={sampleOptions.length < 2}
              name="sample"
              label="Sample Name"
              options={sampleOptions(study, cancer)}
              control={control}
            />
          </Col>
          <Col lg="auto" className="d-flex">
            <Button
              className="mt-auto mb-3"
              disabled={!profile || !sample || !publicSample}
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
      <div id="pcWithinPlot">
        {error && (
          <>
            <hr />
            <div className="p-3">
              <p>An error has occured. Please verify your input.</p>
              <p>{error.data}</p>
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
