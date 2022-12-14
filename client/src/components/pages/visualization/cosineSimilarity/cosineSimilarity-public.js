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
import { defaultMatrix2, defaultMatrix } from '../../../../services/utils';
import { useSeqmatrixOptionsQuery } from '../../../../services/store/rootApi';

export default function CsPublic() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.visualization);
  const mergeForm = (state) =>
    dispatch(
      actions.mergeVisualization({
        cosineSimilarity: { ...store.cosineSimilarity, withinForm: state },
      })
    );

  const { matrixData, matrixList, projectID } = store.main;
  const { withinForm } = store.cosineSimilarity;

  const [params, setParams] = useState('');

  const { data, error, isFetching } = useCosineWithinQuery(params, {
    skip: !params,
  });

  const { data: publicOptions } = useSeqmatrixOptionsQuery();

  const { control, handleSubmit, watch, setValue } = useForm({
    defaultValues: withinForm,
  });

  const { profile, matrix, study, cancer } = watch();

  const profileOptions = matrixData.length
    ? [...new Set(matrixData.map((e) => e.profile))].map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  const matrixOptions = (profile) =>
    matrixData && profile
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

  const studyOptions = publicOptions.length
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
    if (!profile && profileOptions.length) handleProfile(profileOptions[0]);
  }, [profileOptions]);
  useEffect(() => {
    if (!study && studyOptions.length) handleStudy(studyOptions[0]);
  }, [studyOptions]);

  function onSubmit(data) {
    mergeForm(data);
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
      projectID,
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

  const customStyles = {
    control: (base, state) => ({
      ...base,
      background: '#f1e4ef',
      // match with the menu
      borderRadius: state.isFocused ? '3px 3px 0 0' : 3,
      // Overwrittes the different states of border
      borderColor: state.isFocused ? '#f1e4ef' : '#8e4b86',
      // Removes weird border around container
      boxShadow: state.isFocused ? null : null,
      '&:hover': {
        // Overwrittes the different states of border
        borderColor: state.isFocused ? '#8e4b86' : '#f1e4ef',
      },
    }),
    menu: (base) => ({
      ...base,
      // override border radius to match the box
      borderRadius: 0,
      // kill the gap
      marginTop: 0,
    }),
    menuList: (base) => ({
      ...base,
      // kill the white space on first and last option
      padding: 0,
    }),
  };

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
              styles={customStyles}
            />
          </Col>
          <Col lg="auto">
            <Select
              disabled={!profile}
              name="matrix"
              label="Matrix Size"
              options={matrixOptions(profile)}
              control={control}
              styles={customStyles}
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
              styles={customStyles}
            />
          </Col>
          <Col lg="auto">
            <Select
              disabled={!study}
              name="cancer"
              label="Cancer Type"
              options={cancerOptions(study)}
              control={control}
              styles={customStyles}
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
              plotPath={'web/results/' + data.output.plotPath}
              txtPath={`web/results/${data.output.plotPath}`}
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
