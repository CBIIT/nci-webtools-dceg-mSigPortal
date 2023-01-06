import { useState, useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { useForm } from 'react-hook-form';
import { NavHashLink } from 'react-router-hash-link';
import Plotly from '../../../../controls/plotly/plot/plot';
import { LoadingOverlay } from '../../../../controls/loading-overlay/loading-overlay';
import Select from '../../../../controls/select/selectForm';
import { actions as catalogActions } from '../../../../../services/store/catalog';
import { actions as modalActions } from '../../../../../services/store/modal';
import { useSignatureOptionsQuery } from '../../../../../services/store/rootApi';
import { useRsComparisonQuery } from './apiSlice';

const actions = { ...catalogActions, ...modalActions };

export default function RsComparisonPlot() {
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
  } = useRsComparisonQuery(params, { skip: !params });

  const { control, setValue, watch, handleSubmit } = useForm({
    defaultValues: cosineSimilarityStore,
  });

  const { profile, matrix, signatureSet1, signatureSet2 } = watch();

  const supportedProfiles = ['SBS', 'DBS', 'ID'];
  const supportedMatrices = [96, 192, 78, 83];
  const profileOptions = data
    ? [...new Set(data.map((e) => e.profile))]
        .filter((e) => supportedProfiles.includes(e))
        .map((e) => ({
          label: e,
          value: e,
        }))
    : [];

  const matrixOptions = (profile) =>
    data && profile
      ? [
          ...new Set(
            data
              .filter((e) => e.profile == profile.value)
              .map((e) => e.matrix)
              .filter((e) => supportedMatrices.includes(e))
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

  const signatureNameOptions = (profile, matrix, signatureSet) =>
    data && profile
      ? [
          ...new Set(
            data
              .filter(
                (e) =>
                  e.profile == profile.value &&
                  e.matrix == matrix.value &&
                  e.signatureSetName == signatureSet.value
              )
              .map((e) => e.signatureName)
              .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
          ),
        ].map((e) => ({
          label: e,
          value: e,
        }))
      : [];

  // set inital parameters
  useEffect(() => {
    if (!profile && profileOptions.length) handleProfile(profileOptions[0]);
  }, [profile, profileOptions]);

  function handleProfile(profile) {
    const matrices = matrixOptions(profile);

    setValue('profile', profile);
    handleMatrix(profile, matrices[0]);
  }

  function handleMatrix(profile, matrix) {
    const signatureSets = signatureSetOptions(profile, matrix);
    const sigSet1 = signatureSets[0];
    const sigSet2 = signatureSets[1] || signatureSets[0];

    setValue('matrix', matrix);
    handleSignatureSet(profile, matrix, sigSet1, 1);
    handleSignatureSet(profile, matrix, sigSet2, 2);
  }

  function handleSignatureSet(profile, matrix, signatureSet, index) {
    const signatureNames = signatureNameOptions(profile, matrix, signatureSet);

    setValue(`signatureSet${index}`, signatureSet);
    setValue(`signatureName${index}`, signatureNames[0]);
  }

  function onSubmit(data) {
    const params = {
      profile: data.profile.value,
      matrix: data.matrix.value,
      signatureSetName: `${data.signatureSet1.value};${data.signatureSet2.value}`,
      signatureName: `${data.signatureName1.value};${data.signatureName2.value}`,
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
              name="matrix"
              label="Matrix Size"
              options={matrixOptions(profile)}
              disabled={fetchingOptions || fetchingPlot}
              onChange={(e) => handleMatrix(profile, e)}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="signatureSet1"
              label="Reference Signature Set 1"
              options={signatureSetOptions(profile, matrix)}
              disabled={fetchingOptions || fetchingPlot}
              onChange={(e) => handleSignatureSet(profile, matrix, e, 1)}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="signatureName1"
              label="Signature Name 1"
              options={signatureNameOptions(profile, matrix, signatureSet1)}
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
              onChange={(e) => handleSignatureSet(profile, matrix, e, 2)}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="signatureName2"
              label="Signature Name 2"
              options={signatureNameOptions(profile, matrix, signatureSet2)}
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

        {plot && !plotError && (
          <>
            <Plotly
              data={plot.traces}
              layout={plot.layout}
              config={plot.config}
              originalData={plot.original}
              filename="cosine_similarity"
            />
            <div className="p-4">
              <p>
                The plot above shows the mutational profiles of two selected
                signatures, as well as the difference between them. The text at
                the top of the plot indicates the profile similarity calculated
                using Residual Sum of Squares (RSS) and cosine similarity
                methods.
              </p>
              <p>
                Residual Sum of Squares (RSS) measures the discrepancy between
                two mutational profiles. Cosine similarity measures how similar
                two mutational profiles are. For example, two identical
                mutational signatures will have RSS = 0 and Cosine similarity =
                1. For additional information about RSS and cosine similarity,
                click{' '}
                <NavHashLink to="/faq#cosine-similarity">here</NavHashLink>.
              </p>
            </div>
          </>
        )}
      </div>
    </div>
  );
}
