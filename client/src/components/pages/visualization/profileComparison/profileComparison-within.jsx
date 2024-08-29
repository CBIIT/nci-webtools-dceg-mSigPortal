import { useState, useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import Select from '../../../controls/select/selectHookForm';
import { useForm } from 'react-hook-form';
import { NavHashLink } from 'react-router-hash-link';
import { useProfileComparisonWithinQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plotly from '../../../controls/plotly/plot/plot';
import { useSeqmatrixOptionsQuery } from '../../../../services/store/rootApi';
import { defaultProfile2 } from '../../../../services/utils';

export default function PcWithin({ state }) {
  const [params, setParams] = useState(null);
  const { study, cancer, strategy, id, source } = state;

  const { data: options } = useSeqmatrixOptionsQuery(
    {
      ...(source == 'public'
        ? { study: study.value, cancer: cancer.value, strategy: strategy.value }
        : { userId: id }),
    },
    { skip: source == 'user' ? !id : !study }
  );
  const { data, error, isFetching } = useProfileComparisonWithinQuery(params, {
    skip: !params,
  });

  const { control, handleSubmit, watch, setValue } = useForm({
    defaultValues: { profile: '', sample1: '', sample2: '' },
  });

  const { profile, sample1, sample2 } = watch();

  const profileOptions = options
    ? [
        ...new Set(
          options
            .map((e) => e.profile)
            .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
        ),
      ].map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  const sampleOptions = (profile) =>
    options && profile
      ? [
          ...new Set(
            options
              .filter((e) => e.profile === profile.value)
              .map((e) => e.sample)
              .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
          ),
        ].map((e) => ({
          label: e,
          value: e,
        }))
      : [];

  // set inital parameters
  // useEffect(() => {
  //   if (!profile && profileOptions.length) handleProfile(profileOptions[0]);
  // }, [profileOptions]);

  useEffect(() => {
    if (!profile && profileOptions.length)
      handleProfile(defaultProfile2(profileOptions));
  }, [profileOptions]);

  function onSubmit(data) {
    const params = {
      ...(source == 'public' && {
        study: study.value,
        cancer: cancer.value,
        strategy: strategy.value,
      }),

      ...(source == 'user' && { userId: id }),

      profile: data.profile.value,
      matrix:
        data.profile.value === 'SBS'
          ? '96'
          : data.profile.value === 'DBS'
          ? '78'
          : '83',
      sample: data.sample1.value + ';' + data.sample2.value,
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
              disabled={sampleOptions(profile).length < 2}
              name="sample1"
              label="Sample Name 1"
              options={sampleOptions(profile)}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              disabled={sampleOptions(profile).length < 2}
              name="sample2"
              label="Sample Name 2"
              options={sampleOptions(profile)}
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
        {sampleOptions(profile).length < 2 && (
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
