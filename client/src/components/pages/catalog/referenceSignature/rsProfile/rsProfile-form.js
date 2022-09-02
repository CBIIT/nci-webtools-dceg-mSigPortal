import React, { useEffect } from 'react';
import { Row, Col, Button, Form } from 'react-bootstrap';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import { LoadingOverlay } from '../../../../controls/loading-overlay/loading-overlay';
import Plotly from '../../../../controls/plotly/plot/plot';
import SvgContainer from '../../../../controls/svgContainer/svgContainer';
import CustomSelect from '../../../../controls/select/select-old';
import Description from '../../../../controls/description/description';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../../services/store/catalog';
import { actions as modalActions } from '../../../../../services/store/modal';
import { useForm } from 'react-hook-form';
import Select from '../../../../controls/select/selectForm';
import { useSignatureOptionsQuery } from '../../../../../services/store/rootApi';
import { useRsProfileQuery } from './apiSlice';
import {
  defaultProfile,
  defaultMatrix,
  defaultFilter,
} from '../../../../../services/utils';

const actions = { ...catalogActions, ...modalActions };

export default function Profile({ submitR }) {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.catalog);

  const mergeSigMutationalProfiles = (state) =>
    dispatch(actions.mergeCatalog({ sigMutationalProfiles: state }));

  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const { source, matrixList, projectID } = store.main;
  const { plots, debugR, err, loading } = store.sigMutationalProfiles;
  const {
    refSigData,
    sample,
    profile,
    matrix,
    signatureSet,
    strategy,
    signatureName,
  } = store.referenceSignature;

  const { control, setValue, watch } = useForm();

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
  console.log(signatureSourceOptions);

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
                  e.sample == sample.value &&
                  e.profile == profile.value &&
                  supportMatrix[e.profile].includes(e.matrix)
              )
              .map((e) => e.matrix)
              .sort((a, b) => a - b)
          ),
        ].map((e) => ({ label: e, value: e }))
      : [];

  const referenceSignatureSetOption = (source, profile, matrix) =>
    source && profile && matrix
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
        ]
      : [];

  const strategyOptions = (source, profile, matrix, signatureSet) =>
    source && profile && matrix && signatureSet
      ? [
          ...new Set(
            optiondata.map((e) => e.strategy).sort((a, b) => b.localeCompare(a))
          ),
        ].map((e) => ({ label: e, value: e }))
      : [];

  const signatureNameOptions = (
    source,
    profile,
    matrix,
    signatureSet,
    strategy
  ) =>
    source && profile && matrix && signatureSet && strategy
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
    const profile = defaultProfile(profiles);

    mergeSigMutationalProfiles({
      source: source,
      profile: profile,
      matrix: matrix,
      signatureSet: signatureSet,
    });
  }

  function handleProfile(profile) {
    const matrices = matrixOptions(sample, profile);
    const matrix = defaultMatrix(profile, matrices);
    const signatureSet = referenceSignatureSetOption(source, profile, matrix);
    const strategy = strategyOptions(source, profile, matrix, signatureSet);
    const signatureName = signatureNameOptions(
      source,
      profile,
      matrix,
      signatureSet,
      strategy
    );
    mergeSigMutationalProfiles({
      profile: profile,
      matrix: matrix,
      signatureSet: signatureSet,
      strategy: strategy,
      signatureName: signatureName,
    });
  }

  function handleMatrix(matrix) {
    const signatureSet = referenceSignatureSetOption(source, profile, matrix);
    const strategy = strategyOptions(source, profile, matrix, signatureSet);
    const signatureName = signatureNameOptions(
      source,
      profile,
      matrix,
      signatureSet,
      strategy
    );
    mergeSigMutationalProfiles({
      matrix: matrix,
      signatureSet: signatureSet,
      strategy: strategy,
      signatureName: signatureName,
    });
  }

  function handleSet(signatureSet) {
    const strategy = strategyOptions(signatureSet);
    const signatureName = signatureNameOptions(
      source,
      profile,
      matrix,
      signatureSet,
      strategy
    );
    mergeSigMutationalProfiles({
      signatureSet: signatureSet,
      strategy: strategy,
      signatureName: signatureName,
    });
  }

  function handleStrategy(trategy) {
    const signatureName = signatureNameOptions(trategy);

    mergeSigMutationalProfiles({
      strategy: strategy,
      signatureName: signatureName,
    });
  }

  function handleName(signatureName) {
    mergeSigMutationalProfiles({
      signatureName: signatureName,
    });
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
            <CustomSelect
              id="mspSource"
              label="Signature Source"
              value={source}
              options={signatureSourceOptions}
              control={control}
              onChange={handleSource}
            />
          </Col>
          <Col lg="auto">
            <CustomSelect
              id="mspProfileName"
              label="Profile Name"
              value={profile}
              options={profileOptions}
              onChange={handleProfile}
            />
          </Col>
          <Col lg="auto">
            <CustomSelect
              id="mspMatrix"
              label="Matrix"
              value={matrix}
              options={matrixOptions}
              onChange={handleMatrix}
            />
          </Col>
          <Col lg="auto">
            <CustomSelect
              id="mspSet"
              label="Reference Signature Set"
              value={signatureSet}
              options={referenceSignatureSetOption}
              onChange={handleSet}
            />
          </Col>
          <Col lg="auto">
            <CustomSelect
              id="mspStrategy"
              label="Experimental Strategy"
              value={strategy}
              options={strategyOptions}
              onChange={handleStrategy}
            />
          </Col>
          <Col lg="auto">
            <CustomSelect
              id="mspSigName"
              label="Signature Name"
              value={signatureName}
              options={signatureNameOptions}
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
