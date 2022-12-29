import React, { useEffect, useState } from 'react';
import { Row, Col, Button, Form } from 'react-bootstrap';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import { LoadingOverlay } from '../../../../controls/loading-overlay/loading-overlay';
import Plotly from '../../../../controls/plotly/plot/plot';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../../services/store/catalog';
import { actions as modalActions } from '../../../../../services/store/modal';
import { useForm, useFieldArray } from 'react-hook-form';
import Select from '../../../../controls/select/selectForm';
import { useRsProfileOptionsQuery, useRsProfilePlotQuery } from './apiSlice';
import {
  defaultProfile2,
  defaultMatrix2,
  defaultSignatureSet,
  defaultStrategy,
  defaultSignatureName,
} from '../../../../../services/utils';

const actions = { ...catalogActions, ...modalActions };

export default function ProfileFormPlot() {
  const dispatch = useDispatch();
  const { plots: storePlots } = useSelector(
    (state) => state.catalog.rSProfiles
  );

  const mergeState = (state) =>
    dispatch(actions.mergeCatalog({ rSProfiles: state }));

  const [params, setParams] = useState(null);

  const { control, setValue, watch } = useForm({
    defaultValues: { plots: storePlots },
  });

  const {
    fields: plotsFields,
    append: addPlots,
    remove: removePlots,
  } = useFieldArray({
    control,
    name: 'plots',
  });

  //console.log(plotsFields);

  const { plots } = watch();

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
    //source && signatureOptions.length
    source
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
    source && profile && signatureOptions.length
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
    source && profile && matrix && signatureOptions.length
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
            // .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
          ),
        ].map((e) => ({ label: e, value: e }))
      : [];

  const strategyOptions = (source, profile, matrix, signatureSetName) =>
    source && profile && matrix && signatureSetName && signatureOptions.length
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
    source &&
    profile &&
    matrix &&
    signatureSetName &&
    strategy &&
    signatureOptions.length
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

  function handleSource(source, index) {
    const profiles = profileOptions(source);
    const profile = defaultProfile2(profiles);
    const matrices = matrixOptions(source, profile);
    const matrix = defaultMatrix2(profile, matrices);
    const signatureSetNames = referenceSignatureSetOption(
      source,
      profile,
      matrix
    );
    const signatureSetName = defaultSignatureSet(signatureSetNames);
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

    setValue(`plots[${index}].source`, source);
    handleProfile(source, profile, index);
    // setValue(`plots[${index}].profile`, profile);
    // setValue(`plots[${index}].matrix`, matrix);
    // setValue(`plots[${index}].signatureSetName`, signatureSetName);
    // setValue(`plots[${index}].strategy`, strategy);
    // setValue(`plots[${index}].signatureName`, signatureName);

    // mergeState({
    //   source: source,
    //   profile: profile,
    //   matrix: matrix,
    //   signatureSetName: signatureSetName,
    //   strategy: strategy,
    //   signatureName: signatureName,
    // });
  }

  function handleProfile(source, profile, index) {
    const matrices = matrixOptions(source, profile);
    const matrix = defaultMatrix2(profile, matrices);
    const signatureSetNames = referenceSignatureSetOption(
      source,
      profile,
      matrix
    );
    const signatureSetName = defaultSignatureSet(signatureSetNames);
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

    setValue(`plots[${index}].profile`, profile);
    handleMatrix(source, profile, matrix, index);
    //setValue(`plots[${index}].matrix`, matrix);
    // setValue(`plots[${index}].signatureSetName`, signatureSetName);
    // setValue(`plots[${index}].strategy`, strategy);
    // setValue(`plots[${index}].signatureName`, signatureName);
    // mergeState({
    //   profile,
    //   matrix,
    //   signatureSetName,
    //   strategy,
    //   signatureName,
    // });
  }

  function handleMatrix(source, profile, matrix, index) {
    const signatureSetNames = referenceSignatureSetOption(
      source,
      profile,
      matrix
    );
    const signatureSetName = defaultSignatureSet(signatureSetNames);
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

    setValue(`plots[${index}].matrix`, matrix);
    handleSet(source, profile, matrix, signatureSetName, index);
    // setValue(`plots[${index}].signatureSetName`, signatureSetName);
    // setValue(`plots[${index}].strategy`, strategy);
    // setValue(`plots[${index}].signatureName`, signatureName);
    // mergeState({
    //   matrix,
    //   signatureSetName,
    //   strategy,
    //   signatureName,
    // });
  }

  function handleSet(source, profile, matrix, signatureSetName, index) {
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

    setValue(`plots[${index}].signatureSetName`, signatureSetName);
    handleStrategy(source, profile, matrix, signatureSetName, strategy, index);
    // setValue(`plots[${index}].strategy`, strategy);
    // setValue(`plots[${index}].signatureName`, signatureName);

    // mergeState({
    //   signatureSetName,
    //   strategy,
    //   signatureName,
    // });
  }

  function handleStrategy(
    source,
    profile,
    matrix,
    signatureSetName,
    strategy,
    index
  ) {
    const signatureNames = signatureNameOptions(
      source,
      profile,
      matrix,
      signatureSetName,
      strategy
    );
    const signatureName = defaultSignatureName(signatureNames);

    setValue(`plots[${index}].strategy`, strategy);
    handleName(
      source,
      profile,
      matrix,
      signatureSetName,
      strategy,
      signatureName,
      index
    );
    //setValue(`plots[${index}].signatureName`, signatureName);

    // mergeState({
    //   strategy,
    //   signatureName,
    // });
  }

  function handleName(
    source,
    profile,
    matrix,
    signatureSetName,
    strategy,
    signatureName,
    index
  ) {
    setValue(`plots[${index}].signatureName`, signatureName);

    mergeState({
      signatureName,
    });
  }
  // set inital source
  useEffect(() => {
    if (!plots[0].source && signatureSourceOptions.length)
      handleSource(signatureSourceOptions[0]);
  }, [signatureSourceOptions]);

  const validDataArray = plots.filter(
    (e) =>
      e.source?.value &&
      e.profile?.value &&
      e.matrix?.value &&
      e.signatureSetName?.value &&
      e.strategy?.value &&
      e.signatureName?.value
  );
  //console.log(validDataArray);
  //console.log(plots.filter((e) => e.source.value));

  // get data on form change
  useEffect(() => {
    const params = {
      source: plots.filter((e) => e.source.value),
      profile: plots.filter((e) => e.profile.value),
      matrix: plots.filter((e) => e.matrix.value),
      signatureSetName: plots.filter((e) => e.value),
      strategy: plots.filter((e) => e.strategy.value),
      signatureName: plots.filter((e) => e.signatureName.value),
    };
    // const params = plots.filter((e) => ({
    //   source: e.source.value,
    //   profile: e.profile.value,
    //   matrix: e.matrix.value,
    //   signatureSetName: e.signatureSetName.value,
    //   strategy: e.strategy.value,
    //   signatureName: e.signatureName.value,
    // }));
    // console.log(params);
    // console.log(Object.values(params));
    // setParams(params);
    setParams(Object.values(params)[0]);
  }, [plots]);

  const {
    data: plotdata,
    error: plotError,
    isFetching: plotFetching,
  } = useRsProfilePlotQuery(params, {
    skip: !params,
  });
  //console.log(plotdata);
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
      {plotsFields.map((item, index) => (
        <div key={item.id}>
          <Form className="p-3">
            <LoadingOverlay active={plotFetching} />

            <Row className="">
              <Col lg="auto">
                <Select
                  name={`plots[${index}].source`}
                  label="Signature Source"
                  //value={source}
                  options={signatureSourceOptions}
                  control={control}
                  onChange={(e) => handleSource(e, index)}
                  styles={customStyles}
                />
              </Col>
              <Col lg="auto">
                <Select
                  name={`plots[${index}].profile`}
                  label="Profile Name"
                  //value={profile}
                  options={profileOptions(plots[index].source)}
                  control={control}
                  onChange={(e) => handleProfile(plots[index].source, e, index)}
                  styles={customStyles}
                />
              </Col>
              <Col lg="auto">
                <Select
                  name={`plots[${index}] .matrix`}
                  label="Matrix"
                  //value={matrix}
                  options={matrixOptions(
                    plots[index].source,
                    plots[index].profile
                  )}
                  control={control}
                  onChange={(e) =>
                    handleMatrix(
                      plots[index].source,
                      plots[index].profile,
                      e,
                      index
                    )
                  }
                  styles={customStyles}
                />
              </Col>
              <Col lg="auto">
                <Select
                  name={`plots[${index}].signatureSetName`}
                  label="Reference Signature Set"
                  //value={signatureSetName}
                  options={referenceSignatureSetOption(
                    plots[index].source,
                    plots[index].profile,
                    plots[index].matrix
                  )}
                  control={control}
                  onChange={(e) =>
                    handleSet(
                      plots[index].source,
                      plots[index].profile,
                      plots[index].matrix,
                      e,
                      index
                    )
                  }
                  styles={customStyles}
                />
              </Col>
              <Col lg="auto">
                <Select
                  name={`plots[${index}].strategy`}
                  label="Experimental Strategy"
                  //value={strategy}
                  options={strategyOptions(
                    plots[index].source,
                    plots[index].profile,
                    plots[index].matrix,
                    plots[index].signatureSetName
                  )}
                  control={control}
                  onChange={(e) =>
                    handleStrategy(
                      plots[index].source,
                      plots[index].profile,
                      plots[index].matrix,
                      plots[index].signatureSetName,
                      e,
                      index
                    )
                  }
                  styles={customStyles}
                />
              </Col>
              <Col lg="auto">
                <Select
                  name={`plots[${index}].signatureName`}
                  label="Signature Name"
                  //value={signatureName}
                  options={signatureNameOptions(
                    plots[index].source,
                    plots[index].profile,
                    plots[index].matrix,
                    plots[index].signatureSetName,
                    plots[index].strategy
                  )}
                  control={control}
                  onChange={(e) =>
                    handleName(
                      plots[index].source,
                      plots[index].profile,
                      plots[index].matrix,
                      plots[index].signatureSetName,
                      plots[index].strategy,
                      e,
                      index
                    )
                  }
                  styles={customStyles}
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

          <div id={`plots[${index}]`}>
            <div
              style={{ display: plotError || optionError ? 'block' : 'none' }}
            >
              <p>An error has occured. Please verify your input.</p>
            </div>

            {plotdata && (
              <Plotly
                data={plotdata[index].data.traces}
                layout={plotdata[index].data.layout}
                config={plotdata[index].data.config}
                divId="mutationalProfilePlot"
                filename={plots[index].source?.value || 'Mutational Profile'}
              />
            )}
            <Row className="mr-3">
              {index === plotsFields.length - 1 ? (
                <Col className="d-flex justify-content-end">
                  <Button
                    className="ml-auto"
                    variant="link"
                    onClick={() =>
                      addPlots({
                        source: '',
                        profile: '',
                        matrix: '',
                        signatureSetName: '',
                        strategy: '',
                        signatureName: '',
                      })
                    }
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
      ))}
    </div>
  );
}