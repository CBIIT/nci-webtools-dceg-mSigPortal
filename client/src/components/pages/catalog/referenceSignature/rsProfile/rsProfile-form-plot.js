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
import { useSignatureOptionsQuery } from '../../../../../services/store/rootApi';
import { useRsProfileOptionsQuery, useRsProfileDataQuery } from './apiSlice';
import { defaultProfile2, defaultMatrix2 } from '../../../../../services/utils';

const actions = { ...catalogActions, ...modalActions };

export default function ProfileFormPlot({ options, index }) {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.catalog);

  const mergeRsProfiles = (state) =>
    dispatch(actions.mergeCatalog({ rSProfiles: state }));

  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const { matrixList, projectID } = store.main;
  const { plots, debugR, err, loading } = store.rSProfiles;
  const { refSigData, sample } = store.referenceSignature;
  console.log(store);

  const [params, setParams] = useState(null);

  const defaultValues = {
    source: { label: 'Published_signature', value: 'Published_signature' },
    profile: { label: 'SBS', value: 'SBS' },
    matrix: { label: '96', value: '96' },
    signatureSetName: {
      label: 'Other_Published_Signatures_GRCh37_SBS96',
      value: 'Other_Published_Signatures_GRCh37_SBS96',
    },
    strategy: { label: 'WGS', value: 'WGS' },
    signatureName: {
      label: 'SBS_5-FU_organoids_WGS',
      value: 'BS_5-FU_organoids_WGS',
    },
  };
  // const defaultValues = {
  //   source: '',
  //   profile: '',
  //   matrix: '',
  //   signatureSetName: '',
  //   strategy: '',
  //   signatureName: '',
  // };
  const { control, setValue, watch } = useForm();
  //const { control, setValue, watch } = useForm();

  const { source, profile, matrix, signatureSetName, strategy, signatureName } =
    watch();

  console.log(watch());

  const supportMatrix = {
    SBS: [6, 24, 96, 192, 288, 384, 1536],
    DBS: [78, 186],
    ID: [28, 83, 415],
    RS: [32],
    CN: [48],
  };

  const {
    data: optiondata,
    error: optionError,
    isFetching: optionFetching,
  } = useRsProfileOptionsQuery();
  console.log(optiondata);

  const signatureSourceOptions = optiondata
    ? [...new Set(optiondata.map((e) => e.source))]
        .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
        .map((e) => ({
          label: e,
          value: e,
        }))
    : [];

  const profileOptions = (source) =>
    source && optiondata.length
      ? [
          ...new Set(
            optiondata
              .filter((e) => e.source === source.value)
              .map((e) => e.profile)
              .sort((a, b) => a.localeCompare(b))
          ),
        ].map((e) => ({ label: e, value: e }))
      : [];

  const matrixOptions = (source, profile) =>
    source && profile && optiondata.length
      ? [
          ...new Set(
            optiondata
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
    source && profile && matrix && optiondata.length
      ? [
          ...new Set(
            optiondata
              .filter(
                (e) =>
                  e.source === source.value &&
                  e.profile === profile.value &&
                  supportMatrix[e.profile].includes(e.matrix)
              )
              .map((e) => e.signatureSetName)
            // .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
          ),
        ].map((e) => ({ label: e, value: e }))
      : [];

  const strategyOptions = (source, profile, matrix, signatureSetName) =>
    source && profile && matrix && signatureSetName && optiondata.length
      ? [
          ...new Set(
            optiondata
              .filter(
                (e) =>
                  e.source === source.value &&
                  e.profile === profile.value &&
                  supportMatrix[e.profile].includes(e.matrix) &&
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
    source &&
    profile &&
    matrix &&
    signatureSetName &&
    strategy &&
    optiondata.length
      ? [
          ...new Set(
            optiondata
              .filter(
                (e) =>
                  e.source === source.value &&
                  e.profile === profile.value &&
                  supportMatrix[e.profile].includes(e.matrix) &&
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
    const signatureSetName = referenceSignatureSetOption(
      source,
      profile,
      matrix
    );
    const strategy = strategyOptions(source, profile, matrix, signatureSetName);
    const signatureName = signatureNameOptions(
      source,
      profile,
      matrix,
      signatureSetName,
      strategy
    );
    console.log(source);
    console.log(profile);
    console.log(matrices);
    console.log(matrix);
    console.log(signatureSetName);
    console.log(strategy);
    console.log(signatureName);
    setValue('source', source);
    setValue('profile', profile);
    setValue('matrix', matrix);
    setValue('signatureSetName', signatureSetName);
    setValue('strategy', strategy);
    setValue('signatureName', signatureName);
    // mergeRsProfiles({
    //   source,
    //   profile,
    //   matrix,
    //   signatureSetName,
    //   strategy,
    //   signatureName,
    // });
  }

  function handleProfile(profile) {
    const matrices = matrixOptions(source, profile);
    const matrix = defaultMatrix2(profile, matrices);
    const signatureSetName = referenceSignatureSetOption(
      source,
      profile,
      matrix
    );
    const strategy = strategyOptions(source, profile, matrix, signatureSetName);
    const signatureName = signatureNameOptions(
      source,
      profile,
      matrix,
      signatureSetName,
      strategy
    );
    console.log(profile);
    console.log(matrix);
    console.log(matrices);
    console.log(signatureSetName);
    console.log(strategy);
    console.log(signatureName);

    setValue('profile', profile);
    setValue('matrix', matrix);
    setValue('signatureSetName', signatureSetName);
    setValue('strategy', strategy);
    setValue('signatureName', signatureName);
    // mergeRsProfiles({
    //   profile,
    //   matrix,
    //   signatureSetName,
    //   strategy,
    //   signatureName,
    // });
  }

  function handleMatrix(matrix) {
    const signatureSetName = referenceSignatureSetOption(
      source,
      profile,
      matrix
    );
    const strategy = strategyOptions(source, profile, matrix, signatureSetName);
    const signatureName = signatureNameOptions(
      source,
      profile,
      matrix,
      signatureSetName,
      strategy
    );
    console.log(profile);
    console.log(matrix);
    console.log(signatureSetName);
    console.log(strategy);
    console.log(signatureName);

    setValue('matrix', matrix);
    setValue('signatureSetName', signatureSetName);
    setValue('strategy', strategy);
    setValue('signatureName', signatureName);
    // mergeRsProfiles({
    //   matrix,
    //   signatureSetName,
    //   strategy,
    //   signatureName,
    // });
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

    console.log(profile);
    console.log(matrix);
    console.log(signatureSetName);
    console.log(strategy);
    console.log(signatureName);

    setValue('signatureSetName', signatureSetName);
    setValue('strategy', strategy);
    setValue('signatureName', signatureName);

    // mergeRsProfiles({
    //   signatureSetName,
    //   strategy,
    //   signatureName,
    // });
  }

  function handleStrategy(strategy) {
    const signatureName = signatureNameOptions(strategy);
    setValue('strategy', strategy);
    setValue('signatureName', signatureName);
    // mergeRsProfiles({
    //   strategy,
    //   signatureName,
    // });
  }

  function handleName(signatureName) {
    setValue('signatureName', signatureName);
    // mergeRsProfiles({
    //   signatureName,
    // });
  }
  // set inital source
  useEffect(() => {
    if (!source && signatureSourceOptions.length)
      handleSource(signatureSourceOptions[0]);
  }, [signatureSourceOptions]);

  console.log(matrixOptions[0]);
  console.log(signatureSetName);

  // set inital signature set name
  // useEffect(() => {
  //   if (!signatureSetName && referenceSignatureSetOption.length)
  //     handleSet(referenceSignatureSetOption[0]);
  // }, [referenceSignatureSetOption]);
  // get data on form change
  useEffect(() => {
    console.log(source);
    console.log(profile);
    console.log(matrix);
    console.log(signatureSetName);
    console.log(strategy);
    console.log(signatureName);
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
      console.log(params);
      setParams(params);
    }
  }, [source, profile, matrix, signatureSetName, strategy, signatureName]);

  const {
    data: plotdata,
    error: plotError,
    isFetching: plotFetching,
  } = useRsProfileDataQuery(params, {
    skip: !params,
  });
  console.log(plotdata);

  function addPlots() {
    mergeRsProfiles({
      plots: [
        ...plots,
        {
          source: '',
          profile: '',
          matrix: '',
          signatureSetName: '',
          strategy: '',
          signatureName: '',
        },
      ],
    });
  }
  function removePlots(index) {
    //if (plots[index.plotURL]) Object.revokeObjectURL(plots[index].plotURL);
    let newPlots = plots.slice();
    newPlots.splice(index, 1);

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
        {index === 0 ? (
          <Row className="mt-3">
            <Col md="auto" className="d-flex">
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
          </Row>
        ) : (
          <Row className="mt-3">
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
                </span>
              </Button>
            </Col>
          </Row>
        )}
      </Form>

      <div id="plot0">
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
      </div>
      {/* {additionalPlots()} */}
    </div>
  );
}
