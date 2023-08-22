import { useState, useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import Select from '../../../controls/select/selectHookForm';
import { useForm } from 'react-hook-form';
import { NavHashLink } from 'react-router-hash-link';
import {
  useCosineReferenceQuery,
  useCosineSignatureSetsQuery,
} from './apiSlice';
import Description from '../../../controls/description/description';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import SvgContainer from '../../../controls/svgContainer/svgContainer';
import { defaultProfile2, defaultMatrix } from '../../../../services/utils';
import { useSeqmatrixOptionsQuery } from '../../../../services/store/rootApi';
import { useMatrixListQuery } from '../userForm/apiSlice';

export default function CsReference({ state }) {
  const [params, setParams] = useState('');
  const [signatureSetQuery, setSignatureSetQuery] = useState('');
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
  // plot
  const { data, error, isFetching } = useCosineReferenceQuery(params, {
    skip: !params,
  });
  // get signature sets
  const { data: signatureSetOptions, isFetching: fetchingSigSets } =
    useCosineSignatureSetsQuery(signatureSetQuery, {
      skip: !signatureSetQuery,
    });

  // main form
  const { control, handleSubmit, watch, setValue } = useForm({
    defaultValues: { profile: '', signatureSet: '' },
  });
  const { profile, signatureSet } = watch();

  // declare form Options
  const profileOptions = options
    ? [...new Set(options.map((e) => e.profile))].map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  // set initial profile
  useEffect(() => {
    if (!profile && profileOptions.length)
      setValue('profile', defaultProfile2(profileOptions));
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
    const { profile, signatureSet } = data;
    const cacheBust = new Date().getTime();
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
            cacheBust,
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
            cacheBust,
          };
    setParams(params);
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
              plotPath={'api/data/' + data.output.plotPath}
              txtPath={`api/data/${data.output.plotPath}`}
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
