import { useState, useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import Select from '../../../controls/select/selectHookForm';
import { useForm } from 'react-hook-form';
import { NavHashLink } from 'react-router-hash-link';
import { useCosinePublicQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import SvgContainer from '../../../controls/svgContainer/svgContainer';
import Description from '../../../controls/description/description';
import {
  defaultMatrix2,
  defaultMatrix,
  defaultProfile2,
} from '../../../../services/utils';
import { useSeqmatrixOptionsQuery } from '../../../../services/store/rootApi';
import { useMatrixListQuery } from '../userForm/apiSlice';

export default function CsPublic({ state }) {
  const { id } = state;

  const [params, setParams] = useState('');

  const { data, error, isFetching } = useCosinePublicQuery(params, {
    skip: !params,
  });
  const { data: options } = useSeqmatrixOptionsQuery(
    { userId: id },
    { skip: !id }
  );
  const { data: matrixList } = useMatrixListQuery(id, { skip: !id });
  const { data: publicOptions } = useSeqmatrixOptionsQuery({
    columns: ['study', 'cancer'],
  });

  const { control, handleSubmit, watch, setValue } = useForm({
    defaultValues: { profile: '', matrix: '', study: '', cancer: '' },
  });

  const { profile, matrix, study, cancer } = watch();

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

  const studyOptions = publicOptions
    ? [...new Set(publicOptions.map((e) => e.study))].map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  const cancerOptions = (study) =>
    publicOptions && study
      ? [
          ...new Set(
            publicOptions
              .filter((e) => e.study == study.value)
              .map((e) => e.cancer)
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
  useEffect(() => {
    if (!study && studyOptions.length) handleStudy(studyOptions[0]);
  }, [studyOptions]);

  function onSubmit(data) {
    const cacheBust = new Date().getTime();
    const params = {
      fn: 'cosineSimilarityPublic',
      args: {
        profileName: data.profile.value + data.matrix.value,
        study: data.study.value,
        cancerType: data.cancer.value,
        matrixFile: matrixList.filter(
          (e) =>
            e.profile == data.profile.value &&
            e.matrix == defaultMatrix(data.profile.value, ['96', '78', '83'])
        )[0].Path,
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

  function handleStudy(e) {
    const cancers = cancerOptions(e);
    setValue('study', e);
    if (cancers.length) {
      setValue('cancer', cancers[0]);
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
              disabled={!profile}
              name="matrix"
              label="Matrix Size"
              options={matrixOptions(profile)}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              disabled={!studyOptions.length}
              name="study"
              label="Study"
              options={studyOptions}
              onChange={handleStudy}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              disabled={!study}
              name="cancer"
              label="Cancer Type"
              options={cancerOptions(study)}
              control={control}
            />
          </Col>
          <Col lg="auto" className="d-flex">
            <Button
              className="mt-auto mb-3"
              disabled={!profile || !matrix || !study || !cancer}
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
              plotPath={'api/data/' + data.output.plotPath}
              txtPath={`api/data/${data.output.plotPath}`}
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
