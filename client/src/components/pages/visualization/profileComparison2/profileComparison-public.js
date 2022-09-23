import { useState, useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import Select from '../../../controls/select/selectForm';
import { useForm } from 'react-hook-form';
import { NavHashLink } from 'react-router-hash-link';
import { useSelector, useDispatch } from 'react-redux';
import { actions } from '../../../../services/store/visualization';
import { useProfileComparisonWithinQuery } from './apiSlice';
import { useSeqmatrixOptionsQuery } from '../publicForm/apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import SvgContainer from '../../../controls/svgContainer/svgContainer';
import { defaultMatrix } from '../../../../services/utils';

export default function PcWithin() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.visualization);
  const mergeForm = (state) =>
    dispatch(
      actions.mergeVisualization({
        profileComparison: { ...store.profileComparison, withinForm: state },
      })
    );

  const { study, cancer, strategy } = store.publicForm;
  const { source, svgList, matrixList, projectID } = store.main;
  const { withinForm } = store.profileComparison;

  const [params, setParams] = useState(null);

  //   const { data: publicData } = useSeqmatrixOptionsQuery({
  //     skip: source == 'user',
  //   });
  const { data, error, isFetching } = useProfileComparisonWithinQuery(params, {
    skip: !params,
  });

  const getSampleOptions = (profile) =>
    svgList && profile
      ? [
          ...new Set(
            svgList
              .filter((e) => e.profile == profile.value)
              .map((e) => e.sample)
          ),
        ].map((e) => ({
          label: e,
          value: e,
        }))
      : [];

  const { control, handleSubmit, watch, setValue } = useForm({
    defaultValues: withinForm,
  });

  const { profile, sample1, sample2 } = watch();

  const profileOptions = svgList.length
    ? [...new Set(svgList.map((e) => e.profile))].map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  const sampleOptions = getSampleOptions(profile);

  // set inital parameters
  useEffect(() => {
    if (!profile && profileOptions.length) handleProfile(profileOptions[0]);
  }, [profileOptions]);

  function onSubmit(data) {
    mergeForm(data);

    const params =
      source == 'user'
        ? {
            fn: 'profileComparisonWithin',
            args: {
              profileType: data.profile.value,
              sampleName1: data.sample1.value,
              sampleName2: data.sample2.value,
              matrixFile: matrixList.filter(
                (row) =>
                  row.profileType == data.profile.value &&
                  row.matrixSize ==
                    defaultMatrix(data.profile.value, ['96', '78', '83'])
              )[0].Path,
            },
            projectID,
          }
        : {
            fn: 'profileComparisonWithinPublic',
            args: {
              profileType: data.profile.value,
              sampleName1: data.sample1.value,
              sampleName2: data.sample2.value,
              study: study.value,
              cancerType: cancer.value,
              experimentalStrategy: strategy.value,
            },
            projectID,
          };
    setParams(params);
  }

  function handleProfile(e) {
    const samples = getSampleOptions(e);

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
        Input a [Profile Type] and two sample names ([Sample Name 1], [Sample
        Name 2]) to generate the mutational profile of each sample, as well as
        the difference between the two mutational profiles.
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
              disabled={sampleOptions.length < 2}
              name="sample1"
              label="Sample Name 1"
              options={sampleOptions}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              disabled={sampleOptions.length < 2}
              name="sample2"
              label="Sample Name 2"
              options={sampleOptions}
              control={control}
            />
          </Col>
          <Col lg="auto" className="d-flex">
            <Button
              className="mt-auto mb-3"
              disabled={!profile || !sample1 || !sample2}
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
        {data?.output.plotPath && (
          <>
            <hr />
            <SvgContainer
              className="p-3"
              downloadName={data.output.plotPath.split('/').slice(-1)[0]}
              plotPath={'web/results/' + data.output.plotPath}
              txtPath={`web/results/${data.output.plotPath}`}
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
