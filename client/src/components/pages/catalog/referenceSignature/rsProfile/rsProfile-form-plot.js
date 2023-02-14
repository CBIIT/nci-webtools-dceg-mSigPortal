import React, { useEffect, useState } from 'react';
import { Row, Col, Button, Form } from 'react-bootstrap';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import { LoadingOverlay } from '../../../../controls/loading-overlay/loading-overlay';
import Plotly from '../../../../controls/plotly/plot/plot';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../../services/store/catalog';
import { actions as modalActions } from '../../../../../services/store/modal';
import { useForm } from 'react-hook-form';
import Select from '../../../../controls/select/selectForm';
import { useRsProfileOptionsQuery, useRsProfilePlotQuery } from './apiSlice';
import {
  defaultProfile2,
  defaultMatrix2,
  defaultSignatureSet2,
  defaultStrategy,
  defaultSignatureName,
} from '../../../../../services/utils';

const actions = { ...catalogActions, ...modalActions };

export default function ProfileFormPlot({ options, index }) {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.catalog);

  const mergeRsProfiles = (state) =>
    dispatch(actions.mergeCatalog({ rSProfiles: state }));

  const mergeState = (state) => {
    let newPlot = plots.slice();
    newPlot[index] = { ...newPlot[index], ...state };
    mergeRsProfiles({ plots: newPlot });
  };

  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const { matrixList, id } = store.main;
  const { refSigData, sample } = store.referenceSignature;

  const { plots, err, loading } = store.rSProfiles;

  const [params, setParams] = useState(null);

  const defaultValues = {
    source: '',
    profile: '',
    matrix: '',
    signatureSetName: '',
    strategy: '',
    signatureName: '',
  };
  const { control, setValue, watch } = useForm({ defaultValues: plots[index] });
  //const { control, setValue, watch } = useForm({ defaultValues });
  const { source, profile, matrix, signatureSetName, strategy, signatureName } =
    watch();

  const supportMatrix = {
    SBS: [6, 24, 96, 192, 288, 384, 1536],
    DBS: [78, 186],
    ID: [28, 83, 415],
    RS: [32],
    CN: [48],
  };

  const {
    data: signatureOptions,
    error: optionError,
    isFetching: optionFetching,
  } = useRsProfileOptionsQuery();

  const signatureSourceOptions = signatureOptions
    ? [...new Set(signatureOptions.map((e) => e.source))]
        .sort((a, b) => b.localeCompare(a, undefined, { numeric: true }))
        .map((e) => ({
          label: e,
          value: e,
        }))
    : [];

  const profileOptions = (source) =>
    signatureOptions && source
      ? [
          ...new Set(
            signatureOptions
              .filter((e) => e.source === source.value)
              .map((e) => e.profile)
              .sort((a, b) => a.localeCompare(b))
          ),
        ].map((e) => ({ label: e, value: e }))
      : [];

  const matrixOptions = (source, profile) =>
    signatureOptions && source && profile
      ? [
          ...new Set(
            signatureOptions
              .filter(
                (e) =>
                  e.source === source.value &&
                  e.profile === profile.value &&
                  supportMatrix[e.profile].includes(e.matrix)
              )
              .map((e) => e.matrix)
              .sort((a, b) => a - b)
          ),
        ].map((e) => ({ label: e, value: e }))
      : [];

  const referenceSignatureSetOption = (source, profile, matrix) =>
    signatureOptions && source && profile && matrix
      ? [
          ...new Set(
            signatureOptions
              .filter(
                (e) =>
                  e.source === source.value &&
                  e.profile === profile.value &&
                  e.matrix === matrix.value
              )
              .map((e) => e.signatureSetName)
              .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
          ),
        ].map((e) => ({ label: e, value: e }))
      : [];

  const strategyOptions = (source, profile, matrix, signatureSetName) =>
    signatureOptions && source && profile && matrix && signatureSetName
      ? [
          ...new Set(
            signatureOptions
              .filter(
                (e) =>
                  e.source === source.value &&
                  e.profile === profile.value &&
                  e.matrix === matrix.value &&
                  e.signatureSetName === signatureSetName.value
              )
              .map((e) => e.strategy)
              .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
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
    signatureOptions &&
    source &&
    profile &&
    matrix &&
    signatureSetName &&
    strategy
      ? [
          ...new Set(
            signatureOptions
              .filter(
                (e) =>
                  e.source === source.value &&
                  e.profile === profile.value &&
                  e.matrix === matrix.value &&
                  e.signatureSetName === signatureSetName.value &&
                  e.strategy === strategy.value
              )
              .map((e) => e.signatureName)
              .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
          ),
        ].map((e) => ({ label: e, value: e }))
      : [];

  function handleSource(source) {
    const profiles = profileOptions(source);
    const profile = defaultProfile2(profiles);
    const matrices = matrixOptions(source, profile);
    const matrix = defaultMatrix2(profile, matrices);
    const signatureSetNames = referenceSignatureSetOption(
      source,
      profile,
      matrix
    );
    const signatureSetName = defaultSignatureSet2(signatureSetNames);
    const strategies = strategyOptions(
      source,
      profile,
      matrix,
      signatureSetName
    );
    const strategy = defaultStrategy(strategies);
    const signatureNames = signatureNameOptions(
      source,
      profile,
      matrix,
      signatureSetName,
      strategy
    );
    const signatureName = defaultSignatureName(signatureNames);

    setValue('source', source);
    setValue('profile', profile);
    setValue('matrix', matrix);
    setValue('signatureSetName', signatureSetName);
    setValue('strategy', strategy);
    setValue('signatureName', signatureName);
    mergeState({
      source: source,
      profile: profile,
      matrix: matrix,
      signatureSetName: signatureSetName,
      strategy: strategy,
      signatureName: signatureName,
    });
  }

  function handleProfile(profile) {
    const matrices = matrixOptions(source, profile);
    const matrix = defaultMatrix2(profile, matrices);
    const signatureSetNames = referenceSignatureSetOption(
      source,
      profile,
      matrix
    );
    const signatureSetName = defaultSignatureSet2(signatureSetNames);
    const strategies = strategyOptions(
      source,
      profile,
      matrix,
      signatureSetName
    );
    const strategy = defaultStrategy(strategies);
    const signatureNames = signatureNameOptions(
      source,
      profile,
      matrix,
      signatureSetName,
      strategy
    );
    const signatureName = defaultSignatureName(signatureNames);

    setValue('profile', profile);
    setValue('matrix', matrix);
    setValue('signatureSetName', signatureSetName);
    setValue('strategy', strategy);
    setValue('signatureName', signatureName);
    mergeState({
      profile,
      matrix,
      signatureSetName,
      strategy,
      signatureName,
    });
  }

  function handleMatrix(matrix) {
    const signatureSetNames = referenceSignatureSetOption(
      source,
      profile,
      matrix
    );
    const signatureSetName = defaultSignatureSet2(signatureSetNames);
    const strategies = strategyOptions(
      source,
      profile,
      matrix,
      signatureSetName
    );
    const strategy = defaultStrategy(strategies);
    const signatureNames = signatureNameOptions(
      source,
      profile,
      matrix,
      signatureSetName,
      strategy
    );
    const signatureName = defaultSignatureName(signatureNames);

    setValue('matrix', matrix);
    setValue('signatureSetName', signatureSetName);
    setValue('strategy', strategy);
    setValue('signatureName', signatureName);
    mergeState({
      matrix,
      signatureSetName,
      strategy,
      signatureName,
    });
  }

  function handleSet(signatureSetName) {
    const strategies = strategyOptions(
      source,
      profile,
      matrix,
      signatureSetName
    );
    const strategy = defaultStrategy(strategies);
    const signatureNames = signatureNameOptions(
      source,
      profile,
      matrix,
      signatureSetName,
      strategy
    );
    const signatureName = defaultSignatureName(signatureNames);

    setValue('signatureSetName', signatureSetName);
    setValue('strategy', strategy);
    setValue('signatureName', signatureName);

    mergeState({
      signatureSetName,
      strategy,
      signatureName,
    });
  }

  function handleStrategy(strategy) {
    const signatureNames = signatureNameOptions(
      source,
      profile,
      matrix,
      signatureSetName,
      strategy
    );
    const signatureName = defaultSignatureName(signatureNames);
    setValue('strategy', strategy);
    setValue('signatureName', signatureName);

    mergeState({
      strategy,
      signatureName,
    });
  }

  function handleName(signatureName) {
    setValue('signatureName', signatureName);
    mergeState({
      signatureName,
    });
  }
  // set inital source
  useEffect(() => {
    if (!source && signatureSourceOptions.length)
      handleSource(signatureSourceOptions[0]);
  }, [signatureSourceOptions]);

  // get data on form change
  useEffect(() => {
    if (
      source?.value &&
      profile?.value &&
      matrix?.value &&
      signatureSetName?.value &&
      strategy?.value &&
      signatureName?.value
    ) {
      const params = {
        source: source.value,
        profile: profile.value,
        matrix: matrix.value,
        signatureSetName: signatureSetName.value,
        strategy: strategy.value,
        signatureName: signatureName.value,
      };

      setParams(params);
    }
  }, [source, profile, matrix, signatureSetName, strategy, signatureName]);

  const {
    data: plotdata,
    error: plotError,
    isFetching: plotFetching,
  } = useRsProfilePlotQuery(params, {
    skip: !params,
  });
  console.log(plotdata);

  // function addPlots() {
  //   mergeRsProfiles({
  //     plots: [
  //       ...plots,
  //       {
  //         source: '',
  //         profile: '',
  //         matrix: '',
  //         signatureSetName: '',
  //         strategy: '',
  //         signatureName: '',
  //       },
  //     ],
  //   });
  // }
  //console.log(refSigData);
  function addPlots() {
    const signatureSource = {
      label: 'Reference_signatures',
      value: 'Reference_signatures',
    };
    const profiles = profileOptions(signatureSource);

    const profile = defaultProfile2(profiles);

    const matrices = matrixOptions(source, profile);
    const matrix = defaultMatrix2(profile, matrices);
    const signatureSetNames = referenceSignatureSetOption(
      source,
      profile,
      matrix
    );
    const signatureSetName = defaultSignatureSet2(signatureSetNames);
    const strategies = strategyOptions(
      source,
      profile,
      matrix,
      signatureSetName
    );
    const strategy = defaultStrategy(strategies);
    const signatureNames = signatureNameOptions(
      source,
      profile,
      matrix,
      signatureSetName,
      strategy
    );
    const signatureName = defaultSignatureName(signatureNames);

    mergeRsProfiles({
      plots: [
        ...plots,
        {
          source: signatureSource,
          profile: profile,
          matrix: matrix,
          signatureSetName: signatureSetName,
          strategy: strategy,
          signatureName: signatureName,
          index: plots.length,
        },
      ],
    });
  }

  function removePlots(index) {
    //if (plots[index.plotURL]) Object.revokeObjectURL(plots[index].plotURL);

    let newPlots = plots.slice();

    //newPlots.splice(index, 1);

    const removed = newPlots.splice(index, 1);

    mergeRsProfiles({ plots: newPlots });
  }

  return (
    <div>
      <Form className="p-3">
        <LoadingOverlay active={loading} />
        <Row className="">
          <Col lg="auto">
            <Select
              name="source"
              label="Signature Source"
              value={source}
              options={signatureSourceOptions}
              control={control}
              onChange={handleSource}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="profile"
              label="Profile Name"
              value={profile}
              options={profileOptions(source)}
              control={control}
              onChange={handleProfile}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="matrix"
              label="Matrix"
              value={matrix}
              options={matrixOptions(source, profile)}
              control={control}
              onChange={handleMatrix}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="signatureSetName"
              label="Reference Signature Set"
              value={signatureSetName}
              options={referenceSignatureSetOption(source, profile, matrix)}
              control={control}
              onChange={handleSet}
            />
          </Col>
          <Col lg="auto">
            <Select
              name="strategy"
              label="Experimental Strategy"
              value={strategy}
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
              name="signatureName"
              label="Signature Name"
              value={signatureName}
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
          {index != 0 ? (
            <Col md="auto" className="d-flex">
              <Button
                className="ml-auto"
                variant="link"
                onClick={() => removePlots(index)}
                title={'Remove Plot ' + (parseInt(index) + 1)}
                style={{ textDecoration: 'none' }}
              >
                <span className="text-nowrap" title="Remove Plot">
                  <FontAwesomeIcon icon={faMinus} /> Remove Plot
                </span>{' '}
                {/* {parseInt(index) + 1} */}
              </Button>
            </Col>
          ) : (
            <Col md="auto" className="d-flex"></Col>
          )}
        </Row>
      </Form>

      <div id="plot">
        <div style={{ display: err ? 'block' : 'none' }}>
          <p>An error has occured. Please verify your input.</p>
        </div>

        {plotdata && (
          <Plotly
            data={plotdata.traces}
            layout={plotdata.layout}
            config={plotdata.config}
            divId="mutationalProfilePlot"
            filename={source?.value || 'Mutational Profile'}
          />
        )}
        <Row className="mr-3">
          {index === plots.length - 1 ? (
            <Col className="d-flex justify-content-end">
              <Button
                className="ml-auto"
                variant="link"
                onClick={() => addPlots()}
                title="Add Plot"
                style={{ textDecoration: 'none' }}
              >
                <span className="text-nowrap" title="Add Plot">
                  <FontAwesomeIcon icon={faPlus} /> Add Plot
                </span>
              </Button>
            </Col>
          ) : (
            <Col></Col>
          )}
        </Row>
      </div>

      <hr></hr>
      {/* {additionalPlots()} */}
    </div>
  );
}
