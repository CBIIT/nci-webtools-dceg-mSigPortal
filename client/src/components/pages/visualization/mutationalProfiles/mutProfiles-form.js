import React, { useEffect } from 'react';
import { Form, Row, Col } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import Select from '../../../controls/select/selectForm';
import Description from '../../../controls/description/description';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../../services/store/visualization';
import {
  defaultProfile2,
  defaultMatrix2,
  defaultFilter2,
} from '../../../../services/utils';
import { NavHashLink } from 'react-router-hash-link';
const actions = { ...visualizationActions };

export default function TreeLeafForm() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.visualization);

  const mergeState = (state) =>
    dispatch(actions.mergeVisualization({ mutationalProfiles: state }));

  const { matrixData, source } = store.main;
  const { sample, profile, matrix, filter } = store.mutationalProfiles;

  const { control } = useForm();

  const supportMatrix = {
    SBS: [6, 24, 96, 192, 288, 384, 1536],
    DBS: [78, 186],
    ID: [28, 83, 415],
    RS: [32],
    CN: [48],
  };

  // handle external sample change events
  useEffect(() => {
    if (typeof sample === 'string') {
      const sampleOption =
        sampleOptions.find((e) => e.value == sample) || sampleOptions[0];
      handleSample(sampleOption);
    }
  }, [sample]);

  // populate controls
  useEffect(() => {
    if (matrixData.length && !sample) {
      handleSample(sampleOptions[0]);
    }
  }, [matrixData]);

  const sampleOptions = matrixData.length
    ? [...new Set(matrixData.map((d) => d.sample))]
        .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
        .map((e) => ({
          label: e,
          value: e,
        }))
    : [];

  const profileOptions = (sample) =>
    sample && matrixData.length
      ? [
          ...new Set(
            matrixData
              .filter((e) => e.sample == sample.value)
              .map((e) => e.profile)
              .sort((a, b) => b.localeCompare(a))
          ),
        ].map((e) => ({ label: e, value: e }))
      : [];

  const matrixOptions = (sample, profile) =>
    sample && profile && matrixData.length
      ? [
          ...new Set(
            matrixData
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
    sample && profile && matrix && matrixData.length
      ? [
          ...new Set(
            matrixData
              .filter(
                (e) =>
                  e.sample == sample.value &&
                  e.profile == profile.value &&
                  supportMatrix[e.profile].includes(matrix.value)
              )
              .map((e) => e.filter)
              .sort((a, b) => a - b)
          ),
        ].map((e) => ({ label: e ? e : 'N/A', value: e }))
      : [];

  function handleSample(sample) {
    const profiles = profileOptions(sample);
    const profile = defaultProfile2(profiles);
    const matrices = matrixOptions(sample, profile);
    const matrix = defaultMatrix2(profile, matrices);
    const filters = filterOptions(sample, profile, matrix);
    const filter = defaultFilter2(filters);

    mergeState({ sample, profile, matrix, filter });
  }

  function handleProfile(profile) {
    const matrices = matrixOptions(sample, profile);
    const matrix = defaultMatrix2(profile, matrices);
    const filters = filterOptions(sample, profile, matrix);
    const filter = defaultFilter2(filters);

    mergeState({ profile, matrix, filter });
  }

  function handleMatrix(matrix) {
    const filters = filterOptions(sample, profile, matrix);
    const filter = defaultFilter2(filters);

    mergeState({ matrix, filter });
  }

  function handleFilter(filter) {
    mergeState({ filter: filter });
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
              value={sample}
              options={sampleOptions}
              control={control}
              onChange={handleSample}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="profile"
              label="Profile Type"
              value={profile}
              options={profileOptions(sample)}
              control={control}
              onChange={handleProfile}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="matrix"
              label="Matrix Size"
              value={matrix}
              options={matrixOptions(sample, profile)}
              control={control}
              onChange={handleMatrix}
            />
          </Col>
          {source == 'user' && (
            <Col lg="auto">
              <Select
                name="filter"
                label="Filter"
                value={filter}
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
