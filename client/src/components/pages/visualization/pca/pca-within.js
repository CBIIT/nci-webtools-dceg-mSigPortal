import { useState, useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import Select from '../../../controls/select/selectForm';
import { useForm } from 'react-hook-form';
import { useSelector, useDispatch } from 'react-redux';
import { actions } from '../../../../services/store/visualization';
import { usePcaWithinQuery, usePcaSignatureSetsQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import SvgContainer from '../../../controls/svgContainer/svgContainer';
import { defaultMatrix } from '../../../../services/utils';

export default function CsWithin() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.visualization);
  const mergeForm = (state) =>
    dispatch(
      actions.mergeVisualization({
        pca: { ...store.pca, withinForm: state },
      })
    );

  const { study, cancer, strategy } = store.publicForm;
  const { source, matrixData, matrixList, id } = store.main;
  const { withinForm } = store.pca;

  const [calculationQuery, setCalculationQuery] = useState('');
  const [signatureSetQuery, setSignatureSetQuery] = useState('');

  // get signature sets
  const { data: signatureSetOptions, isFetching: fetchingSigSets } =
    usePcaSignatureSetsQuery(signatureSetQuery, {
      skip: !signatureSetQuery,
    });

  const { data, error, isFetching } = usePcaWithinQuery(calculationQuery, {
    skip: !calculationQuery,
  });

  const { control, handleSubmit, watch, setValue } = useForm({
    defaultValues: withinForm,
  });

  const { profile, signatureSet } = watch();

  const profileOptions = matrixData.length
    ? [...new Set(matrixData.map((e) => e.profile))].map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  // set inital parameters
  useEffect(() => {
    if (!profile && profileOptions.length)
      setValue('profile', profileOptions[0]);
  }, [profileOptions]);
  // set intital signature set
  useEffect(() => {
    if (signatureSetOptions) setValue('signatureSet', signatureSetOptions[0]);
  }, [signatureSetOptions]);

  // get signature sets when profile is selected
  useEffect(() => {
    if (profile) {
      setSignatureSetQuery({
        profile: profile.value,
        matrix: defaultMatrix(profile.value, ['96', '78', '83']),
      });
    }
  }, [profile]);

  function onSubmit(data) {
    mergeForm(data);

    const params =
      source == 'user'
        ? {
            fn: 'pca',
            args: {
              profileType: data.profile.value,
              signatureSet: data.signatureSet.value,
              matrixFile: matrixList.filter(
                (row) =>
                  row.profile == data.profile.value &&
                  row.matrix ==
                    defaultMatrix(data.profile.value, ['96', '78', '83'])
              )[0].Path,
            },
            id,
          }
        : {
            fn: 'pcaPublic',
            args: {
              profileType: data.profile.value,
              signatureSet: data.signatureSet.value,
              study: study.value,
              cancerType: cancer.value,
              experimentalStrategy: strategy.value,
            },
            id,
          };
    setCalculationQuery(params);
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
              control={control}
            />
          </Col>
          <Col lg="auto">
            <Select
              disabled={fetchingSigSets}
              name="signatureSet"
              label="Reference Signature Set"
              options={signatureSetOptions}
              control={control}
            />
          </Col>
          <Col lg="auto" className="d-flex">
            <Button
              className="mt-auto mb-3"
              disabled={!profile || !signatureSet || fetchingSigSets}
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

      {data?.output.heatmap && (
        <div id="heatmapPlot">
          <hr />
          <SvgContainer
            className="p-3"
            downloadName={data.output.heatmap.split('/').slice(-1)[0]}
            plotPath={'web/data/' + data.output.heatmap}
            txtPath={`web/data/${data.output.heatmapData}`}
          />
          <p className="p-3">
            The heatmap shows cosine similarity between each principal component
            and each mutational signature in the selected reference signature
            set. Brighter colors denote higher levels of cosine similarity
            between the principal component and the mutational signature. Red
            dots in some of the boxes indicate a cosine similarity less than 0,
            denoting a negative correlation (e.g., a box with cosine similarity
            value of 0.8 and marked with a red dot indicates a true cosine
            similarity of -0.8).
          </p>
        </div>
      )}
    </div>
  );
}
