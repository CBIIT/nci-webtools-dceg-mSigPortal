import { useState, useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import Select from '../../../controls/select/selectForm';
import { useForm } from 'react-hook-form';
import { useSelector, useDispatch } from 'react-redux';
import { actions } from '../../../../services/store/visualization';
import { usePcaPublicQuery } from './apiSlice';
import { useSeqmatrixOptionsQuery } from '../../../../services/store/rootApi';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import SvgContainer from '../../../controls/svgContainer/svgContainer';

export default function PcaPublic() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.visualization);
  const mergeForm = (state) =>
    dispatch(
      actions.mergeVisualization({
        pca: { ...store.pca, publicForm: state },
      })
    );

  const { matrixData, matrixList, id } = store.main;
  const { publicForm } = store.pca;

  const [params, setParams] = useState(null);

  const { data, error, isFetching } = usePcaPublicQuery(params, {
    skip: !params,
  });
  const { data: publicOptions, isFetching: fetchingPublicOptions } =
    useSeqmatrixOptionsQuery();

  // define form
  const { control, handleSubmit, watch, setValue } = useForm({
    defaultValues: publicForm,
  });
  const { profile, matrix, study, cancer } = watch();

  // define options
  const profileOptions = matrixData.length
    ? [...new Set(matrixData.map((e) => e.profile))].map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  const matrixOptions = (profile) =>
    profile && matrixData.length
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

  const studyOptions = publicOptions
    ? [...new Set(publicOptions.map((e) => e.study))].sort().map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  const cancerOptions = (study) => {
    if (publicOptions && study?.value) {
      return [
        ...[
          ...new Set(
            publicOptions
              .filter((e) => e.study == study.value)
              .map((e) => e.cancer)
          ),
        ]
          .sort()
          .map((e) => ({
            label: e,
            value: e,
          })),
      ];
    } else {
      return [];
    }
  };

  // set inital parameters
  useEffect(() => {
    if (!profile && profileOptions.length) handleProfile(profileOptions[0]);
  }, [profileOptions]);
  useEffect(() => {
    if (!study && studyOptions.length) handleStudy(studyOptions[0]);
  }, [studyOptions]);

  // define form handlers
  function handleProfile(profile) {
    const matrices = matrixOptions(profile);
    setValue('profile', profile);
    setValue('matrix', matrices[0]);
  }

  function handleStudy(study) {
    const cancers = cancerOptions(study);
    setValue('study', study);
    setValue('cancer', cancers[0]);
  }

  function onSubmit(data) {
    const params = {
      fn: 'pcaWithPublic',
      args: {
        profileName: data.profile.value + data.matrix.value,
        study: data.study.value,
        cancerType: data.cancer.value,
        matrixFile: matrixList.filter(
          (row) =>
            row.profile == data.profile.value && row.matrix == data.matrix.value
        )[0].Path,
      },
      id,
    };
    mergeForm(data);
    setParams(params);
  }

  return (
    <div>
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
              disabled={isFetching}
              name="matrix"
              label="Matrix Size"
              options={matrixOptions(profile)}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              disabled={fetchingPublicOptions}
              name="study"
              label="Study"
              options={studyOptions}
              onChange={handleStudy}
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              disabled={fetchingPublicOptions}
              name="cancer"
              label="Cancer Type or Group"
              options={cancerOptions(study)}
              control={control}
            />
          </Col>
          <Col lg="auto" className="d-flex">
            <Button
              className="mt-auto mb-3"
              disabled={!profile || !matrix || !study || !cancer || isFetching}
              variant="primary"
              type="submit"
            >
              Calculate
            </Button>
          </Col>
        </Row>
      </Form>
      {error && (
        <>
          <hr />
          <p className="p-3">An error has occured. Please verify your input.</p>
        </>
      )}

      {data?.output.pca1 && (
        <div id="pca1Plot">
          <hr />
          <SvgContainer
            className="p-3"
            downloadName={data.output.pca1.split('/').slice(-1)[0]}
            plotPath={'web/data/' + data.output.pca1}
          />
          <p className="p-3">
            The bar plot illustrates each of the principal components on the
            x-axis and the percentage of variation that each component explains
            on the y-axis.
          </p>
        </div>
      )}

      {data?.output.pca2 && (
        <div id="data.output.pca2Plot">
          <hr />
          <SvgContainer
            className="p-3"
            downloadName={data.output.pca2.split('/').slice(-1)[0]}
            plotPath={'web/data/' + data.output.pca2}
            txtPath={`web/data/${data.output.pca2Data}`}
          />
          <p className="p-3">
            The individual PCA plot based on the top two principal components
            helps to explain a majority of the variation in selected or input
            data. Each dot on the plot is a sample. The legend on the right
            denotes the percent contribution (contrib) of each sample to the
            principal components on the graph (Dim 1 and Dim 2).
          </p>
        </div>
      )}

      {data?.output.pca3 && (
        <div id="pca3Plot">
          <hr />
          <SvgContainer
            className="p-3"
            downloadName={data.output.pca3.split('/').slice(-1)[0]}
            plotPath={'web/data/' + data.output.pca3}
            txtPath={`web/data/${data.output.pca3Data}`}
          />
          <p className="p-3">
            The variable PCA plot based on the top two principal components
            helps to explain a majority of the variation in the data. The legend
            on the right denotes the percent contribution (contrib) of each
            mutation type to the principal components on the graph (Dim 1 and
            Dim 2).
          </p>
        </div>
      )}
    </div>
  );
}
