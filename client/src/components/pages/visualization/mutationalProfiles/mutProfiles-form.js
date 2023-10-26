import React, { useEffect } from 'react';
import { Form, Row, Col } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import { NavHashLink } from 'react-router-hash-link';
import Select from '../../../controls/select/selectHookForm';
import Description from '../../../controls/description/description';
import {
  defaultProfile2,
  defaultMatrix2,
  defaultFilter2,
} from '../../../../services/utils';
import { useSeqmatrixOptionsQuery } from '../../../../services/store/rootApi';

export default function MutationalProfilesForm({ state, form, mergeForm }) {
  const {
    study,
    cancer,
    strategy,
    id,
    source,
    mutationalProfiles: storeState,
  } = state;
  const { sample: externalSample } = storeState;

  const { data: options } = useSeqmatrixOptionsQuery(
    {
      ...(source == 'public'
        ? { study: study.value, cancer: cancer.value, strategy: strategy.value }
        : { userId: id }),
    },
    { skip: source == 'user' ? !id : !study }
  );
  const { control, setValue, watch } = useForm({ defaultValues: form });
  const { sample, profile, matrix, filter } = watch();

  const supportMatrix = {
    SBS: [6, 24, 96, 192, 288, 384, 1536],
    DBS: [78, 186],
    ID: [28, 83, 415],
    RS: [32],
    CN: [48],
  };

  // handle external sample change events
  useEffect(() => {
    if (typeof externalSample === 'string') {
      const sampleOption =
        sampleOptions.find((e) => e.value == externalSample) ||
        sampleOptions[0];
      handleSample(sampleOption);
    }
  }, [externalSample]);

  // populate controls
  useEffect(() => {
    if (options && options.length && !sample) {
      handleSample(sampleOptions[0]);
    }
  }, [options]);

  // update form
  useEffect(() => {
    if (sample || profile || matrix || filter)
      mergeForm({ sample, profile, matrix, filter });
  }, [sample, profile, matrix, filter]);

  const sampleOptions = options
    ? [...new Set(options.map((d) => d.sample))]
        .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
        .map((e) => ({
          label: e,
          value: e,
        }))
    : [];

  const profileOptions = (sample) =>
    sample && options && options.length
      ? [
          ...new Set(
            options
              .filter((e) => e.sample == sample.value)
              .map((e) => e.profile)
              .sort((a, b) => a.localeCompare(b))
          ),
        ].map((e) => ({ label: e, value: e }))
      : [];

  const matrixOptions = (sample, profile) =>
    sample && profile && options && options.length
      ? [
          ...new Set(
            options
              .filter(
                (e) =>
                  e.sample == sample.value &&
                  e.profile == profile.value &&
                  supportMatrix[e.profile].includes(e.matrix)
              )
              .map((e) => e.matrix)
              .sort((a, b) => a - b)
          ),
        ].map((e) => ({ label: e, value: e }))
      : [];

  const filterOptions = (sample, profile, matrix) =>
    sample && profile && matrix && options && options.length
      ? [
          ...new Set(
            options
              .filter(
                (e) =>
                  e.sample == sample.value &&
                  e.profile == profile.value &&
                  e.matrix == matrix.value
              )
              .map((e) => e.filter)
              .sort((a, b) => a - b)
          ),
        ].map((e) => ({ label: e ? e : 'N/A', value: e }))
      : [];

  function handleSample(sample) {
    const profiles = profileOptions(sample);
    const profile = defaultProfile2(profiles);

    setValue('sample', sample);
    handleProfile(sample, profile);
  }

  function handleProfile(sample, profile) {
    if (profile) {
      const matrices = matrixOptions(sample, profile);
      const matrix = defaultMatrix2(profile, matrices);

      setValue('profile', profile);
      handleMatrix(sample, profile, matrix);
    }
  }

  function handleMatrix(sample, profile, matrix) {
    const filters = filterOptions(sample, profile, matrix);
    const filter = defaultFilter2(filters);

    setValue('matrix', matrix);
    handleFilter(filter);
  }

  function handleFilter(filter) {
    setValue('filter', filter);
  }

  return (
    <div>
      <div className="p-3">
        <Description
          className="m-0"
          less={
            <span>
              Below you can visualize different mutational profiles for a given
              sample.
            </span>
          }
          more={
            <span>
              Use the dropdown arrow to select the [Sample Name], [Profile Type]
              and [Matrix Size]. The [Filter] option is only available if [Split
              Mutations According to Filter] is selected while analyzing user
              data. For additional information on [Profile Type] and [Matrix
              Size], click <NavHashLink to="/faq#sbs">here</NavHashLink>.
            </span>
          }
        />
      </div>
      <hr />
      <Form className="p-3">
        <Row>
          <Col lg="auto">
            <Select
              name="sample"
              label="Sample Name"
              options={sampleOptions}
              control={control}
              onChange={handleSample}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="profile"
              label="Profile Type"
              options={profileOptions(sample)}
              control={control}
              onChange={(e) => handleProfile(sample, e)}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="matrix"
              label="Matrix Size"
              options={matrixOptions(sample, profile)}
              control={control}
              onChange={(e) => handleMatrix(sample, profile, e)}
            />
          </Col>
          {source == 'user' && (
            <Col lg="auto">
              <Select
                name="filter"
                label="Filter"
                disabled={!filterOptions(sample, profile, matrix).length}
                options={filterOptions(sample, profile, matrix)}
                control={control}
                onChange={handleFilter}
              />
            </Col>
          )}
        </Row>
      </Form>
    </div>
  );
}
