import { Button, Form, Row, Col } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { useForm, Controller } from 'react-hook-form';
import Select from '../../../controls/select/selectForm';
import { NavHashLink } from 'react-router-hash-link';
import Description from '../../../controls/description/description';
import { actions as visualizationActions } from '../../../../services/store/visualization';
import {
  defaultProfile2,
  defaultMatrix2,
  defaultFilter2,
} from '../../../../services/utils';
const actions = { ...visualizationActions };

export default function ProfileComparisonForm() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.visualization);

  const mergeState = (state) =>
    dispatch(actions.mergeVisualization({ profileComparison: state }));

  const { matrixData, source } = store.main;
  const { withinForm } = store.profileComparison;

  console.log(store.main);

  console.log(matrixData);

  const { control, handleSubmit, watch, setValue } = useForm({
    defaultValues: withinForm,
  });
  const { profile, sample1, sample2 } = watch();
  const supportMatrix = {
    SBS: [96],
    DBS: [78],
    ID: [83],
  };
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

  async function onSubmit(data) {
    mergeState(data);
  }

  const sampleOptions = matrixData.length
    ? [...new Set(matrixData.map((d) => d.sample))]
        .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
        .map((e) => ({
          label: e,
          value: e,
        }))
    : [];

  const matrixOptions = (sample1, sample2, profile) =>
    sample1 && sample2 && profile && matrixData.length
      ? [
          ...new Set(
            matrixData
              .filter(
                (e) =>
                  e.sample1 &&
                  e.sample2 &&
                  e.profile == profile.value &&
                  supportMatrix[e.profile].includes(e.matrix)
              )
              .map((e) => e.matrix)
              .sort((a, b) => a - b)
          ),
        ].map((e) => ({ label: e, value: e }))
      : [];

  const filterOptions = (samples, profile, matrix) =>
    samples && profile && matrix && matrixData.length
      ? [
          ...new Set(
            matrixData
              .filter(
                (e) => e.sample && e.Profile == profile.value + matrix.value
              )
              .map((e) => e.Filter)
              .sort((a, b) => a - b)
          ),
        ].map((e) => ({ label: e, value: e }))
      : [];

  function handleProfile(profile) {
    const matrices = matrixOptions(sample1, sample2, profile);
    const matrix = defaultMatrix2(profile, matrices);
    const filters = filterOptions(sample1, sample2, profile, matrix);
    const filter = defaultFilter2(filters);

    mergeState({ profile, sample1, sample2, filter });
  }

  return (
    <div>
      <div className="p-3">
        <div>
          <p className="m-0">
            Input a [Profile Type] and two sample names ([Sample Name 1],
            [Sample Name 2]) to generate the mutational profile of each sample,
            as well as the difference between the two mutational profiles.
          </p>
        </div>

        <hr />

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
      </div>
      <hr />
    </div>
  );
}
