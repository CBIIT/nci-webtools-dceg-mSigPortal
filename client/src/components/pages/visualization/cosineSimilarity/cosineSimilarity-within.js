import { useState, useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import Select from '../../../controls/select/selectHookForm';
import { useForm } from 'react-hook-form';
import { NavHashLink } from 'react-router-hash-link';
import { useCosineWithinQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import SvgContainer from '../../../controls/svgContainer/svgContainer';
import Description from '../../../controls/description/description';
import { defaultMatrix2, defaultProfile2 } from '../../../../services/utils';
import { useSeqmatrixOptionsQuery } from '../../../../services/store/rootApi';
import { useMatrixListQuery } from '../userForm/apiSlice';

export default function CsWithin({ state }) {
  const [params, setParams] = useState('');
  const { study, cancer, strategy, source, id } = state;

  const { data: options } = useSeqmatrixOptionsQuery(
    {
      ...(source == 'public'
        ? { study: study.value, cancer: cancer.value, strategy: strategy.value }
        : { userId: id }),
    },
    { skip: source == 'user' ? !id : !study }
  );
  const { data: matrixList } = useMatrixListQuery(id, { skip: !id });
  const { data, error, isFetching } = useCosineWithinQuery(params, {
    skip: !params,
  });

  const { control, handleSubmit, watch, setValue } = useForm({
    defaultValues: { profile: '', matrix: '' },
  });
  const { profile, matrix } = watch();

  const profileOptions = options
    ? [...new Set(options.map((e) => e.profile))].map((e) => ({
        label: e,
        value: e,
      }))
    : [];
  const matrixOptions = (profile) =>
    options && profile
      ? [
          ...new Set(
            options
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

  // set inital parameters
  useEffect(() => {
    if (!profile && profileOptions.length)
      handleProfile(defaultProfile2(profileOptions));
  }, [profileOptions]);

  function onSubmit(data) {
    const cacheBust = new Date().getTime();
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
            cacheBust,
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
            cacheBust,
          };
    setParams(params);
  }

  function handleProfile(e) {
    const matrices = matrixOptions(e);
    setValue('profile', e);
    if (matrices.length) {
      setValue('matrix', defaultMatrix2(e, matrices));
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
              disabled={!matrixOptions(profile).length}
              name="matrix"
              label="Matrix Size"
              options={matrixOptions(profile)}
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
