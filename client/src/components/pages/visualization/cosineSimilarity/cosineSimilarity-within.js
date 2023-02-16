import { useState, useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import Select from '../../../controls/select/selectForm';
import { useForm } from 'react-hook-form';
import { NavHashLink } from 'react-router-hash-link';
import { useSelector, useDispatch } from 'react-redux';
import { actions } from '../../../../services/store/visualization';
import { useCosineWithinQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import SvgContainer from '../../../controls/svgContainer/svgContainer';
import Description from '../../../controls/description/description';
import { defaultMatrix2 } from '../../../../services/utils';

export default function CsWithin() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.visualization);
  const mergeForm = (state) =>
    dispatch(
      actions.mergeVisualization({
        cosineSimilarity: { ...store.cosineSimilarity, withinForm: state },
      })
    );

  const { study, cancer, strategy } = store.publicForm;
  const { source, matrixData, matrixList, id } = store.main;
  const { withinForm } = store.cosineSimilarity;

  const [params, setParams] = useState('');

  const { data, error, isFetching } = useCosineWithinQuery(params, {
    skip: !params,
  });

  function getMatrixOptions(profile) {
    return matrixData && profile
      ? [
          ...new Set(
            matrixData
              .filter((e) => e.profile == profile.value)
              .map((e) => e.matrix)
          ),
        ]
          .sort((a, b) => a - b)
          .map((e) => ({
            label: e,
            value: e,
          }))
      : [];
  }

  const { control, handleSubmit, watch, setValue } = useForm({
    defaultValues: withinForm,
  });

  const { profile, matrix } = watch();

  const profileOptions = matrixData.length
    ? [...new Set(matrixData.map((e) => e.profile))].map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  const matrixOptions = getMatrixOptions(profile);

  // set inital parameters
  useEffect(() => {
    if (!profile && profileOptions.length) handleProfile(profileOptions[0]);
  }, [profileOptions]);

  function onSubmit(data) {
    mergeForm(data);
    const params =
      source == 'user'
        ? {
            fn: 'cosineSimilarityWithin',
            args: {
              matrixFile: matrixList.filter(
                (e) =>
                  e.profile == data.profile.value &&
                  e.matrix == data.matrix.value
              )[0].Path,
            },
            id,
          }
        : {
            fn: 'cosineSimilarityWithinPublic',
            args: {
              profileType: data.profile.value,
              matrixSize: data.matrix.value,
              study: study.value,
              cancerType: cancer.value,
              experimentalStrategy: strategy.value,
            },
            id,
          };
    setParams(params);
  }

  function handleProfile(e) {
    const matrixOptions = getMatrixOptions(e);

    setValue('profile', e);
    if (matrixOptions.length) {
      setValue('matrix', defaultMatrix2(e, matrixOptions));
    }
  }

  return (
    <div>
      <div className="p-3">
        <Description
          less="Cosine similarity is a measure of the similarity of two
                      matrices, which can be helpful to compare two mutational
                      profiles or signatures."
          more={
            <span>
              Below you can explore cosine similarity between sample profiles
              (CS Between Samples), cosine similarity between sample profiles
              and reference signatures (CS to Reference Signatures), or, if
              using your own data, cosine similarity between profiles from your
              input data and profiles from public data (CS to Public Data).
              Simply use the dropdown menus to select a [Profile Type], [Matrix
              Size], or [Reference Signature Set]. Click here to learn more
              about cosine similarity. Click{' '}
              <NavHashLink to="/faq#cosine-similarity">here</NavHashLink> to
              learn more about cosine similarity.
            </span>
          }
        />
      </div>

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
              disabled={!matrixOptions.length}
              name="matrix"
              label="Matrix Size"
              options={matrixOptions}
              control={control}
            />
          </Col>
          <Col lg="auto" className="d-flex">
            <Button
              className="mt-auto mb-3"
              disabled={!profile || !matrix}
              variant="primary"
              type="submit"
            >
              Calculate
            </Button>
          </Col>
        </Row>
      </Form>
      <div id="csWithinPlot">
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
              plotPath={'web/data/' + data.output.plotPath}
              txtPath={`web/data/${data.output.plotPath}`}
            />
            <div className="p-3">
              <p>
                The heatmap shows pairwise cosine similarity between samples
                from the selected profile type. On the x-axis and y-axis are the
                sample names. This analysis will help to highlight the samples
                within the same cluster that may have similar mutational
                signatures.
              </p>
            </div>
          </>
        )}
      </div>
    </div>
  );
}
