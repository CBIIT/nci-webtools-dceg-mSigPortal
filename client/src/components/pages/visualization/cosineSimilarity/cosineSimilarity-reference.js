import { useState, useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import Select from '../../../controls/select/selectForm';
import { useForm } from 'react-hook-form';
import { NavHashLink } from 'react-router-hash-link';
import { useSelector, useDispatch } from 'react-redux';
import { actions } from '../../../../services/store/visualization';
import {
  useCosineReferenceQuery,
  useCosineSignatureSetsQuery,
} from './apiSlice';
import Description from '../../../controls/description/description';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import SvgContainer from '../../../controls/svgContainer/svgContainer';
import { defaultMatrix } from '../../../../services/utils';

export default function CsReference() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.visualization);
  const mergeForm = (state) =>
    dispatch(
      actions.mergeVisualization({
        cosineSimilarity: { ...store.cosineSimilarity, referenceForm: state },
      })
    );

  const { study, cancer, strategy } = store.publicForm;
  const { source, matrixData, matrixList, id } = store.main;
  const { referenceForm } = store.cosineSimilarity;

  // main form
  const {
    control,
    handleSubmit,
    watch,
    setValue,
    resetField,
    formState: { errors: formErrors },
  } = useForm({
    defaultValues: referenceForm,
  });
  const { profile, signatureSet } = watch();

  const [calculationQuery, setCalculationQuery] = useState('');
  const [signatureSetQuery, setSignatureSetQuery] = useState('');

  // get signature sets
  const { data: signatureSetOptions, isFetching: fetchingSigSets } =
    useCosineSignatureSetsQuery(signatureSetQuery, {
      skip: !signatureSetQuery,
    });

  //   calculate plot
  const { data, error, isFetching } = useCosineReferenceQuery(
    calculationQuery,
    { skip: !calculationQuery }
  );

  // declare form Options
  const profileOptions = matrixData.length
    ? [...new Set(matrixData.map((e) => e.profile))].map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  // set inital profile
  useEffect(() => {
    if (!profile && profileOptions.length)
      setValue('profile', profileOptions[0]);
  }, [profileOptions]);
  // set intital signature set
  // useEffect(() => {
  //   if (signatureSetOptions) setValue('signatureSet', signatureSetOptions[0]);
  // }, [signatureSetOptions]);

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

    const { profile, sample, signatureSet, compare } = data;
    const params =
      source == 'user'
        ? {
            fn: 'cosineSimilarityRefSig',
            args: {
              profileType: profile.value,
              signatureSet: signatureSet.value,
              matrixFile: matrixList.filter(
                (e) =>
                  e.profile == profile.value &&
                  e.matrix == defaultMatrix(profile.value, ['96', '78', '83'])
              )[0].Path,
            },
            id,
          }
        : {
            fn: 'cosineSimilarityRefSigPublic',
            args: {
              profileType: profile.value,
              signatureSet: signatureSet.value,
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
      <div className="p-3 m-0">
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
              disabled={!profile || !signatureSet}
              variant="primary"
              type="submit"
            >
              Calculate
            </Button>
          </Col>
        </Row>
      </Form>
      <div id="csReferencePlot">
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
                The following heatmap shows pairwise cosine similarity between
                the mutational profiles of given samples and the selected
                reference signature set. On the x-axis and y-axis are the
                reference signature names and the sample names, respectively.
                This analysis will identify potential dominant mutational
                signatures in selected samples.
              </p>
            </div>
          </>
        )}
      </div>
    </div>
  );
}
