import React, { useEffect } from 'react';
import { Row, Col, Button, Form } from 'react-bootstrap';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import { LoadingOverlay } from '../../../../controls/loading-overlay/loading-overlay';
import Plotly from '../../../../controls/plotly/plot/plot';
import SvgContainer from '../../../../controls/svgContainer/svgContainer';
import Description from '../../../../controls/description/description';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../../services/store/catalog';
import { actions as modalActions } from '../../../../../services/store/modal';
import { useForm } from 'react-hook-form';
import Select from '../../../../controls/select/selectForm';
import { useSignatureOptionsQuery } from '../../../../../services/store/rootApi';
import { useRsProfileQuery } from './apiSlice';
import { defaultProfile2, defaultMatrix2 } from '../../../../../services/utils';

const actions = { ...catalogActions, ...modalActions };

export default function Profile({ submitR }) {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.catalog);

  const mergeSigMutationalProfiles = (state) =>
    dispatch(actions.mergeCatalog({ sigMutationalProfiles: state }));

  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const { matrixList, projectID } = store.main;
  const { plots, debugR, err, loading } = store.sigMutationalProfiles;
  const { refSigData, sample, signatureName } = store.referenceSignature;

  const defaultValues = {
    source: '',
    profile: '',
    matrix: '',
    signatureSetName: '',
    strategy: '',
    signatureName: '',
  };

  const { control, setValue, watch } = useForm({ defaultValues });

  const { source, profile, matrix, signatureSetName, strategy } = watch();

  console.log(watch());

  const supportMatrix = {
    SBS: [6, 24, 96, 192, 288, 384, 1536],
    DBS: [78, 186],
    ID: [28, 83, 415],
  };

  const {
    data: optiondata,
    error: optionError,
    isFetching: optionFetching,
  } = useSignatureOptionsQuery();
  console.log(optiondata);

  useEffect(() => {
    if (optiondata) {
      handleSource(signatureSourceOptions[0]);
    }
  }, []);

  const signatureSourceOptions = optiondata
    ? [...new Set(optiondata.map((e) => e.source))].sort().map((e) => ({
        label: e,
        value: e,
      }))
    : [];

  const profileOptions = (source) =>
    source
      ? [
          ...new Set(
            optiondata
              .filter((e) => e.source === source.value)
              .map((e) => e.profile)
              .sort((a, b) => b.localeCompare(a))
          ),
        ].map((e) => ({ label: e, value: e }))
      : [];

  const matrixOptions = (source, profile) =>
    source && profile
      ? [
          ...new Set(
            optiondata
              .filter(
                (e) =>
                  e.source == source.value &&
                  e.profile == profile.value &&
                  supportMatrix[e.profile].includes(e.matrix)
              )
              .map((e) => e.matrix)
              .sort((a, b) => a - b)
          ),
        ].map((e) => ({ label: e, value: e }))
      : [];

  const referenceSignatureSetOption = (source, profile, matrix) =>
    //console.log(source, profile, matrix);

    source && profile && matrix
      ? [
          ...new Set(
            optiondata
              .filter(
                (e) =>
                  e.source == source.value &&
                  e.profile == profile.value &&
                  supportMatrix[e.profile].includes(matrix.value)
              )
              .map((e) => e.matrix)
              .sort((a, b) => a - b)
          ),
        ]
      : [];

  const strategyOptions = (source, profile, matrix, signatureSetName) =>
    source && profile && matrix && signatureSetName
      ? [
          ...new Set(
            optiondata
              .filter(
                (e) =>
                  e.source == source.value &&
                  e.profile == profile.value &&
                  supportMatrix[e.profile].includes(matrix.value) &&
                  e.signatureName == source.value
              )
              .map((e) => e.strategy)
              .sort((a, b) => b.localeCompare(a))
          ),
        ].map((e) => ({ label: e, value: e }))
      : [];

  const signatureNameOptions = (
    source,
    profile,
    matrix,
    signatureSetName,
    strategy
  ) =>
    source && profile && matrix && signatureSetName && strategy
      ? [
          ...new Set(
            optiondata
              .map((e) => e.signatureName)
              .sort((a, b) => b.localeCompare(a))
          ),
        ].map((e) => ({ label: e, value: e }))
      : [];

  function handleSource(source) {
    const profiles = profileOptions(source);
    const profile = defaultProfile2(profiles);
    const matrices = matrixOptions(sample, profile);
    const matrix = defaultMatrix2(profile, matrices);
    const signatureSetName = referenceSignatureSetOption(
      source,
      profile,
      matrix
    );
    console.log(signatureSetName);
    const strategy = strategyOptions(source, profile, matrix, signatureSetName);
    const signatureName = signatureNameOptions(
      source,
      profile,
      matrix,
      signatureSetName,
      strategy
    );

    setValue('source', source);
    setValue('profile', profile);
    setValue('matrix', matrix);
    setValue('signatureSetName', signatureSetName);
    setValue('strategy', strategy);
    setValue('signatureName', signatureName);
  }

  function handleProfile(profile) {
    const matrices = matrixOptions(sample, profile);
    const matrix = defaultMatrix2(profile, matrices);
    const signatureSetName = referenceSignatureSetOption(
      source,
      profile,
      matrix
    );
    console.log(signatureSetName);
    const strategy = strategyOptions(source, profile, matrix, signatureSetName);
    const signatureName = signatureNameOptions(
      source,
      profile,
      matrix,
      signatureSetName,
      strategy
    );

    setValue('profile', profile);
    setValue('matrix', matrix);
    setValue('signatureSetName', signatureSetName);
    setValue('strategy', strategy);
    setValue('signatureName', signatureName);
  }

  function handleMatrix(matrix) {
    console.log(matrix);
    const signatureSetName = referenceSignatureSetOption(
      source,
      profile,
      matrix
    );
    console.log(signatureSetName);
    const strategy = strategyOptions(source, profile, matrix, signatureSetName);
    const signatureName = signatureNameOptions(
      source,
      profile,
      matrix,
      signatureSetName,
      strategy
    );

    setValue('matrix', matrix);
    setValue('signatureSetName', signatureSetName);
    setValue('strategy', strategy);
    setValue('signatureName', signatureName);
  }

  function handleSet(signatureSetName) {
    const strategy = strategyOptions(source, profile, matrix, signatureSetName);
    const signatureName = signatureNameOptions(
      source,
      profile,
      matrix,
      signatureSetName,
      strategy
    );

    setValue('signatureSetName', signatureSetName);
    setValue('strategy', strategy);
    setValue('signatureName', signatureName);
  }

  function handleStrategy(trategy) {
    const signatureName = signatureNameOptions(trategy);
    setValue('strategy', strategy);
    setValue('signatureName', signatureName);
  }

  function handleName(signatureName) {
    mergeSigMutationalProfiles({
      signatureName: signatureName,
    });

    setValue('signatureName', signatureName);
  }

  const {
    data: plotdata,
    error: plotError,
    isFetching: plotFetching,
  } = useRsProfileQuery();
  console.log(plotdata);

  return (
    <div>
      <Description
        className="p-3 m-0"
        less="Enter any [Signature Source], [Profile Name], [Reference Signature Set], [Experimental Strategy], and [Signature Name] below to visualize the mutational signature profile."
        more="Click ‘+ Add Plot’ to load one or more mutational signature profiles at the same time."
      />
      <hr />
      <Form className="p-3">
        <LoadingOverlay active={loading} />
        <Row className="">
          <Col lg="auto">
            <Select
              name="source"
              label="Signature Source"
              //value={source}
              options={signatureSourceOptions}
              control={control}
              onChange={handleSource}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="profile"
              label="Profile Name"
              //value={profile}
              options={profileOptions(source)}
              control={control}
              onChange={handleProfile}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="matrix"
              label="Matrix"
              //value={matrix}
              options={matrixOptions(source, profile)}
              control={control}
              onChange={handleMatrix}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="signatureSetName"
              label="Reference Signature Set"
              //value={signatureSet}
              options={referenceSignatureSetOption(source, profile, matrix)}
              control={control}
              onChange={handleSet}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="strategy"
              label="Experimental Strategy"
              //value={strategy}
              options={strategyOptions(
                source,
                profile,
                matrix,
                signatureSetName
              )}
              control={control}
              onChange={handleStrategy}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="signature"
              label="Signature Name"
              //value={signatureName}
              options={signatureNameOptions(
                source,
                profile,
                matrix,
                signatureSetName,
                strategy
              )}
              control={control}
              onChange={handleName}
            />
          </Col>
        </Row>
        {/* <AdditionalControls /> */}
        <Row className="mt-3">
          <Col md="auto" className="d-flex">
            <Button
              className="ml-auto"
              variant="link"
              //onClick={() => addPlots()}
              title="Add Plot"
              style={{ textDecoration: 'none' }}
            >
              <span className="text-nowrap" title="Add Plot">
                <FontAwesomeIcon icon={faPlus} /> Add Plot
              </span>
            </Button>
          </Col>
        </Row>
      </Form>

      <div id="plot0">
        <div style={{ display: err ? 'block' : 'none' }}>
          <p>An error has occured. Please verify your input.</p>
        </div>
        {/* {plots[0].plotURL && (
          <>
            <hr />
            <SvgContainer
              className="p-3"
              //title={plotTitle(plots[0])}
              downloadName={plots[0].plotPath.split('/').slice(-1)[0]}
              plotPath={plots[0].plotURL}
              height="500px"
            />
            
          </>
        )} */}
        {optiondata && (
          <Plotly
            data={optiondata.traces}
            layout={optiondata.layout}
            config={optiondata.config}
            divId="mutationalProfilePlot"
            filename={sample?.value || 'Mutational Profile'}
          />
        )}
      </div>
      {/* {additionalPlots()} */}
    </div>
  );
}
