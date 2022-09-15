import React, { useEffect, useState } from 'react';
import { Row, Col, Button, Form } from 'react-bootstrap';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import { LoadingOverlay } from '../../../../controls/loading-overlay/loading-overlay';
import Plotly from '../../../../controls/plotly/plot/plot';
import Description from '../../../../controls/description/description';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../../services/store/catalog';
import { actions as modalActions } from '../../../../../services/store/modal';
import { useForm } from 'react-hook-form';
import Select from '../../../../controls/select/selectForm';
import { useSignatureOptionsQuery } from '../../../../../services/store/rootApi';
import { useRsProfileOptionsQuery, useRsProfileDataQuery } from './apiSlice';
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
  const { refSigData, sample } = store.referenceSignature;

  const [params, setParams] = useState(null);

  // const defaultValues = {
  //   source: { label: 'Published_signature', value: 'Published_signature' },
  //   profile: { label: 'SBS', value: 'SBS' },
  //   matrix: { label: '96', value: '96' },
  //   signatureSetName: {
  //     label: 'Other_Published_Signatures_GRCh37_SBS96',
  //     value: 'Other_Published_Signatures_GRCh37_SBS96',
  //   },
  //   strategy: { label: 'WGS', value: 'WGS' },
  //   signatureName: { label: '', value: '' },
  // };
  const defaultValues = {
    source: '',
    profile: '',
    matrix: '',
    signatureSetName: '',
    strategy: '',
    signatureName: '',
  };
  const { control, setValue, watch } = useForm({ defaultValues });

  const { source, profile, matrix, signatureSetName, strategy, signatureName } =
    watch();

  console.log(watch());

  const supportMatrix = {
    SBS: [6, 24, 96, 192, 288, 384, 1536],
    DBS: [78, 186],
    ID: [28, 83, 415],
    RS: [32],
  };

  const {
    data: optiondata,
    error: optionError,
    isFetching: optionFetching,
  } = useRsProfileOptionsQuery();
  console.log(optiondata);

  useEffect(() => {
    if (optiondata) {
      handleSource(signatureSourceOptions[0]);
    }
  }, []);

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
                (e) => e.source === source.value && e.profile === profile.value
              )
              .map((e) => e.matrix)
              .sort((a, b) => a - b)
          ),
        ].map((e) => ({ label: e, value: e }))
      : [];

  const referenceSignatureSetOption = (source, profile, matrix) =>
    //console.log(source, profile, matrix);

    source && profile && matrix && optiondata.length
      ? [
          ...new Set(
            optiondata
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
    source && profile && matrix && signatureSetName && optiondata.length
      ? [
          ...new Set(
            optiondata
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
    const matrices = matrixOptions(sample, profile);
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

    setValue('matrix', matrix);
    setValue('signatureSetName', signatureSetName[0]);
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

  function handleStrategy(strategy) {
    const signatureName = signatureNameOptions(strategy);
    setValue('strategy', strategy);
    setValue('signatureName', signatureName);
  }

  function handleName(signatureName) {
    setValue('signatureName', signatureName);
  }

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

  const AdditionalControls = () => {
    let controls = [];
    for (let index in plots) {
      controls.push(
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
              //value={signatureSetName}
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
              name="signatureName"
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
      );
    }
    return controls.slice(1);
  };

  const additionalPlots = () => {
    let display = [];
    for (let index in plots) {
      display.push(
        <div id={'plot' + index} key={'plot' + index}>
          <div style={{ display: err ? 'block' : 'none' }}>
            <p>An error has occured. Please verify your input.</p>
          </div>
          {plots[index].plotURL && (
            <div style={{ display: plots[index].plotURL ? 'block' : 'none' }}>
              <hr />
              <span className="font-weight-bold p-3">
                Plot {parseInt(index) + 1}
              </span>
              <Plotly
                data={plotdata.traces}
                layout={plotdata.layout}
                config={plotdata.config}
                divId="mutationalProfilePlot"
                filename={sample?.value || 'Mutational Profile'}
              />
            </div>
          )}
        </div>
      );
    }
    return display.slice(1);
  };

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
              //value={signatureSetName}
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
              name="signatureName"
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
        <AdditionalControls />
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
        {plotdata && (
          <Plotly
            data={plotdata.traces}
            layout={plotdata.layout}
            config={plotdata.config}
            divId="mutationalProfilePlot"
            filename={sample?.value || 'Mutational Profile'}
          />
        )}
      </div>
      {additionalPlots()}
    </div>
  );
}
