import { useState, useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { useForm } from 'react-hook-form';
import Plotly from '../../../../controls/plotly/plot/plot';
import { LoadingOverlay } from '../../../../controls/loading-overlay/loading-overlay';
import Select from '../../../../controls/select/selectForm';
import { actions as catalogActions } from '../../../../../services/store/catalog';
import { actions as modalActions } from '../../../../../services/store/modal';
import { useSignatureOptionsQuery } from '../../../../../services/store/rootApi';
import { useCosineSimilarityQuery } from './apiSlice';

const actions = { ...catalogActions, ...modalActions };

export default function CosineSimilarityPlot() {
  const dispatch = useDispatch();
  const mergeState = (state) =>
    dispatch(actions.mergeCatalog({ cosineSimilarity: state }));

  const cosineSimilarityStore = useSelector(
    (state) => state.catalog.cosineSimilarity
  );

  const [params, setParams] = useState(false);

  // query parameter options
  const {
    data,
    error,
    isFetching: fetchingOptions,
  } = useSignatureOptionsQuery();
  // query data and generate plot after submit
  const {
    data: plot,
    error: plotError,
    isFetching: fetchingPlot,
  } = useCosineSimilarityQuery(params, { skip: !params });

  const {
    register,
    control,
    setValue,
    watch,
    handleSubmit,
    formState: { errors },
  } = useForm({ defaultValues: cosineSimilarityStore });

  const { profile, matrix } = watch();

  // set inital parameters
  useEffect(() => {
    if (!profile && data) handleProfile(profileOptions[0]);
  }, [profile, data]);
  // const supportedProfiles = ['SBS', 'DBS', 'ID'];
  // const supportedMatrices = [96, 192, 78, 73];

  const profileOptions = data
    ? [...new Set(data.map((e) => e.profile))].map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  const matrixOptions = (profile) =>
    data && profile
      ? [
          ...new Set(
            data.filter((e) => e.profile == profile.value).map((e) => e.matrix)
          ),
        ].map((e) => ({
          label: e,
          value: e,
        }))
      : [];

  const signatureSetOptions = (profile, matrix) =>
    data && profile
      ? [
          ...new Set(
            data
              .filter(
                (e) => e.profile == profile.value && e.matrix == matrix.value
              )
              .map((e) => e.signatureSetName)
              .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
          ),
        ].map((e) => ({
          label: e,
          value: e,
        }))
      : [];

  function handleProfile(profile) {
    const matrices = matrixOptions(profile);

    setValue('profile', profile);
    handleMatrix(profile, matrices[0]);
  }

  function handleMatrix(profile, matrix) {
    const signatureSets = signatureSetOptions(profile, matrix);

    setValue('matrix', matrix);
    setValue('signatureSet1', signatureSets[0]);
    setValue('signatureSet2', signatureSets[1] || signatureSets[0]);
  }

  function onSubmit(data) {
    const params = {
      profile: data.profile.value,
      matrix: data.matrix.value,
      signatureSetName: `${data.signatureSet1.value};${data.signatureSet2.value}`,
    };
    setParams(params);
    mergeState(data);
  }

  return (
    <div>
      <Form className="p-3" onSubmit={handleSubmit(onSubmit)}>
        <LoadingOverlay active={fetchingOptions || fetchingPlot} />
        <Row className="">
          <Col lg="auto">
            <Select
              name="profile"
              label="Profile Name"
              options={profileOptions}
              onChange={handleProfile}
              disabled={fetchingOptions || fetchingPlot}
              control={control}
              
            />
          </Col>
          <Col lg="auto">
            <Select
              control={control}
              name="matrix"
              label="Matrix Size"
              options={matrixOptions(profile)}
              disabled={fetchingOptions || fetchingPlot}
              onChange={(e) => handleMatrix(profile, e)}
              
            />
          </Col>
          <Col lg="auto">
            <Select
              name="signatureSet1"
              label="Reference Signature Set 1"
              options={signatureSetOptions(profile, matrix)}
              disabled={fetchingOptions || fetchingPlot}
              control={control}
              
            />
          </Col>
          <Col lg="auto">
            <Select
              name="signatureSet2"
              label="Reference Signature Set 2"
              options={signatureSetOptions(profile, matrix)}
              disabled={fetchingOptions || fetchingPlot}
              control={control}
              
            />
          </Col>
          <Col lg="auto" className="d-flex justify-content-end">
            <Button
              className="mt-auto mb-3"
              variant="primary"
              type="submit"
              disabled={fetchingOptions || fetchingPlot}
            >
              Calculate
            </Button>
          </Col>
        </Row>
      </Form>
      <hr />
      <div>
        <div
          className="text-center"
          style={{ display: plotError ? 'block' : 'none' }}
        >
          <p>An error has occured. Please verify your input.</p>
        </div>

        {plot && (
          <>
            <Plotly
              data={plot.traces}
              layout={plot.layout}
              config={plot.config}
              originalData={plot.original}
              filename="cosine_similarity"
            />
            <p className="p-3">
              The heatmap above shows the cosine similarities between two
              mutational signature sets given a profile type. The text on the
              bottom and left show the signature names of the two selected
              reference signature sets.
            </p>
          </>
        )}
      </div>
    </div>
  );
}
