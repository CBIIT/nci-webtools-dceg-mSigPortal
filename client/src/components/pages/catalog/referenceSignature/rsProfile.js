import React, { useEffect } from 'react';
import { Row, Col, Button, Form } from 'react-bootstrap';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import SvgContainer from '../../../controls/svgContainer/svgContainer';

import CustomSelect from '../../../controls/select/select-old';
import Description from '../../../controls/description/description';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../services/store/catalog';
import { actions as modalActions } from '../../../../services/store/modal';

const actions = { ...catalogActions, ...modalActions };

export default function Profile({ submitR }) {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.catalog);
  const mergeSigMutationalProfiles = (state) =>
    dispatch(actions.mergeCatalog({ sigMutationalProfiles: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const { plots, debugR, err, loading } = store.sigMutationalProfiles;
  const { refSigData } = store.referenceSignature;

  // get inital plot
  useEffect(() => {
    async function getInitalPlot() {
      const plot = plots[0];
      const path = buildPlotPath(
        plot.profileName,
        plot.rsSet,
        plot.signatureName
      );
      const url = await fetchPlot(path, 0);
      mergeSigMutationalProfiles({
        plots: [{ ...plots[0], plotPath: path, plotURL: url }],
      });
    }

    if (plots.length == 1 && plots[0].profileName && !plots[0].plotURL) {
      getInitalPlot();
    }
  }, [plots]);

  // dataFolder/all_svg/rsSet/profileName+signatureName
  function buildPlotPath(profileName, rsSet, signatureName) {
    const profile = profileName.match(/[a-z]+|\d+/gi).join('_');
    const set = rsSet.replace(/\s/g, '_');
    // s3 key
    return `msigportal/Database/Signature/all_svg/${set}/${profile}_plots_mSigPortal_${signatureName}.svg`;
  }

  async function fetchPlot(path, index = null) {
    try {
      const svg = await (
        await fetch(`web/getImageS3`, {
          method: 'POST',
          headers: {
            Accept: 'image/svg',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({ path: path }),
        })
      ).blob();

      if (index && plots[index].plotURL)
        URL.revokeObjectURL(plots[index].plotURL);
      return URL.createObjectURL(svg);
    } catch (err) {
      return err;
    }
  }

  async function handleSource(source, index) {
    const filteredData = refSigData.filter((row) => row.Source == source);
    const profileNameOptions = [
      ...new Set(filteredData.map((row) => row.Profile)),
    ];
    const profileName = profileNameOptions[0];
    const rsSetOptions = [
      ...new Set(
        filteredData
          .filter((row) => row.Source == source && row.Profile == profileName)
          .map((row) => row.Signature_set_name)
      ),
    ];
    const rsSet = rsSetOptions[0];
    const strategyOptions = [
      ...new Set(
        filteredData
          .filter(
            (row) =>
              row.Source == source &&
              row.Profile == profileName &&
              row.Signature_set_name == rsSet
          )
          .map((row) => row.Dataset)
      ),
    ];
    const strategy = strategyOptions[0];
    const signatureNameOptions = [
      ...new Set(
        filteredData
          .filter(
            (row) =>
              row.Source == source &&
              row.Profile == profileName &&
              row.Signature_set_name == rsSet &&
              row.Dataset == strategy
          )
          .map((row) => row.Signature_name)
      ),
    ];

    const plotPath = buildPlotPath(profileName, rsSet, signatureNameOptions[0]);
    const plotURL = await fetchPlot(plotPath, index);

    let newPlots = plots.slice();
    newPlots[index] = {
      ...newPlots[index],
      plotPath: plotPath,
      plotURL: plotURL,
      signatureSource: source,
      profileName: profileName,
      profileNameOptions: profileNameOptions,
      rsSet: rsSet,
      rsSetOptions: rsSetOptions,
      strategy: strategy,
      strategyOptions: strategyOptions,
      signatureName: signatureNameOptions[0],
      signatureNameOptions: signatureNameOptions,
    };

    mergeSigMutationalProfiles({
      plots: newPlots,
    });
  }

  async function handleProfile(profile, index) {
    const filteredData = refSigData.filter(
      (row) =>
        row.Source == plots[index].signatureSource && row.Profile == profile
    );
    const rsSetOptions = [
      ...new Set(filteredData.map((row) => row.Signature_set_name)),
    ];
    const rsSet = rsSetOptions[0];
    const strategyOptions = [
      ...new Set(
        filteredData
          .filter(
            (row) =>
              row.Source == plots[index].signatureSource &&
              row.Profile == profile &&
              row.Signature_set_name == rsSet
          )
          .map((row) => row.Dataset)
      ),
    ];
    const strategy = strategyOptions[0];
    const signatureNameOptions = [
      ...new Set(
        filteredData
          .filter(
            (row) =>
              row.Source == plots[index].signatureSource &&
              row.Profile == profile &&
              row.Signature_set_name == rsSet &&
              row.Dataset == strategy
          )
          .map((row) => row.Signature_name)
      ),
    ];

    const plotPath = buildPlotPath(profile, rsSet, signatureNameOptions[0]);

    const plotURL = await fetchPlot(plotPath, index);

    let newPlots = plots.slice();
    newPlots[index] = {
      ...newPlots[index],
      plotPath: plotPath,
      plotURL: plotURL,
      profileName: profile,
      rsSet: rsSet,
      rsSetOptions: rsSetOptions,
      strategy: strategy,
      strategyOptions: strategyOptions,
      signatureName: signatureNameOptions[0],
      signatureNameOptions: signatureNameOptions,
    };

    mergeSigMutationalProfiles({
      plots: newPlots,
    });
  }

  async function handleSet(set, index) {
    let filteredData = refSigData.filter(
      (row) =>
        row.Source == plots[index].signatureSource &&
        row.Profile == plots[index].profileName &&
        row.Signature_set_name == set
    );
    const strategyOptions = [
      ...new Set(filteredData.map((row) => row.Dataset)),
    ];
    const strategy = strategyOptions[0];
    const signatureNameOptions = [
      ...new Set(
        filteredData
          .filter(
            (row) =>
              row.Source == plots[index].signatureSource &&
              row.Profile == plots[index].profileName &&
              row.Signature_set_name == set &&
              row.Dataset == strategy
          )
          .map((row) => row.Signature_name)
      ),
    ];

    const plotPath = buildPlotPath(
      plots[index].profileName,
      set,
      signatureNameOptions[0]
    );

    const plotURL = await fetchPlot(plotPath, index);

    let newPlots = plots.slice();
    newPlots[index] = {
      ...newPlots[index],
      plotPath: plotPath,
      plotURL: plotURL,
      rsSet: set,
      strategy: strategy,
      strategyOptions: strategyOptions,
      signatureName: signatureNameOptions[0],
      signatureNameOptions: signatureNameOptions,
    };

    mergeSigMutationalProfiles({
      plots: newPlots,
    });
  }

  async function handleStrategy(strategy, index) {
    let filteredData = refSigData.filter(
      (row) =>
        row.Source == plots[index].signatureSource &&
        row.Profile == plots[index].profileName &&
        row.Signature_set_name == plots[index].rsSet &&
        row.Dataset == strategy
    );

    const signatureNameOptions = [
      ...new Set(filteredData.map((row) => row.Signature_name)),
    ];

    const plotPath = buildPlotPath(
      plots[index].profileName,
      plots[index].rsSet,
      signatureNameOptions[0]
    );

    const plotURL = await fetchPlot(plotPath, index);

    let newPlots = plots.slice();
    newPlots[index] = {
      ...newPlots[index],
      plotPath: plotPath,
      plotURL: plotURL,
      strategy: strategy,
      signatureName: signatureNameOptions[0],
      signatureNameOptions: signatureNameOptions,
    };

    mergeSigMutationalProfiles({
      plots: newPlots,
    });
  }

  async function handleName(name, index) {
    const plotPath = buildPlotPath(
      plots[index].profileName,
      plots[index].rsSet,
      name
    );

    const plotURL = await fetchPlot(plotPath, index);

    let newPlots = plots.slice();
    newPlots[index] = {
      ...newPlots[index],
      plotPath: plotPath,
      plotURL: plotURL,
      signatureName: name,
    };

    mergeSigMutationalProfiles({
      plots: newPlots,
    });
  }

  async function addPlots() {
    const signatureSourceOptions = [
      ...new Set(refSigData.map((row) => row.Source)),
    ];
    const signatureSource = signatureSourceOptions[0];
    const profileNameOptions = [
      ...new Set(
        refSigData
          .filter((row) => row.Source == signatureSource)
          .map((row) => row.Profile)
      ),
    ];
    const profileName = profileNameOptions[0];
    const rsSetOptions = [
      ...new Set(
        refSigData
          .filter(
            (row) => row.Source == signatureSource && row.Profile == profileName
          )
          .map((row) => row.Signature_set_name)
      ),
    ];
    const rsSet = rsSetOptions[0];

    const strategyOptions = [
      ...new Set(
        refSigData
          .filter(
            (row) =>
              row.Source == signatureSource &&
              row.Profile == profileName &&
              row.Signature_set_name == rsSet
          )
          .map((row) => row.Dataset)
      ),
    ];
    const strategy = strategyOptions[0];
    const signatureNameOptions = [
      ...new Set(
        refSigData
          .filter(
            (row) =>
              row.Source == signatureSource &&
              row.Profile == profileName &&
              row.Signature_set_name == rsSet &&
              row.Dataset == strategy
          )
          .map((row) => row.Signature_name)
      ),
    ];

    const plotPath = buildPlotPath(profileName, rsSet, signatureNameOptions[0]);

    const plotURL = await fetchPlot(plotPath, null);

    mergeSigMutationalProfiles({
      plots: [
        ...plots,
        {
          plotPath: plotPath,
          plotURL: plotURL,
          signatureSource: signatureSource,
          signatureSourceOptions: signatureSourceOptions,
          profileName: profileName,
          profileNameOptions: profileNameOptions,
          rsSet: rsSet,
          rsSetOptions: rsSetOptions,
          strategy: strategy,
          strategyOptions: strategyOptions,
          signatureName: signatureNameOptions[0],
          signatureNameOptions: signatureNameOptions,
        },
      ],
    });
  }

  function removePlots(index) {
    if (plots[index.plotURL]) Object.revokeObjectURL(plots[index].plotURL);
    let newPlots = plots.slice();
    newPlots.splice(index, 1);

    mergeSigMutationalProfiles({ plots: newPlots });
  }

  const AdditionalControls = () => {
    let controls = [];
    for (let index in plots) {
      controls.push(
        <Row className="mt-3" key={'control' + index}>
          <Col lg="auto">
            <CustomSelect
              id={`mspSource${index}`}
              label="Signature Source"
              class="msigportal-bg"
              value={plots[index].signatureSource}
              options={plots[index].signatureSourceOptions}
              onChange={(source) => handleSource(source, index)}
              className="msigportal-bg"
            />
          </Col>
          <Col lg="auto">
            <CustomSelect
              id={'mspProfileName' + index}
              label="Profile Name"
              value={plots[index].profileName}
              options={plots[index].profileNameOptions}
              onChange={(profile) => handleProfile(profile, index)}
            />
          </Col>
          <Col lg="auto">
            <CustomSelect
              id={'mspSet' + index}
              label="Reference Signature Set"
              value={plots[index].rsSet}
              options={plots[index].rsSetOptions}
              onChange={(profile) => handleSet(profile, index)}
            />
          </Col>
          <Col lg="auto">
            <CustomSelect
              id={'mspStrategy' + index}
              label="Experimental Strategy"
              value={plots[index].strategy}
              options={plots[index].strategyOptions}
              onChange={(strategy) => handleStrategy(strategy, index)}
            />
          </Col>
          <Col lg="auto">
            <CustomSelect
              id={'mspSigName' + index}
              label="Signature Name"
              value={plots[index].signatureName}
              options={plots[index].signatureNameOptions}
              onChange={(name) => handleName(name, index)}
            />
          </Col>
          <Col lg="auto" className="d-flex">
            <Button
              className="my-auto text-danger"
              variant="link"
              onClick={() => removePlots(index)}
              title={'Remove Plot ' + (parseInt(index) + 1)}
            >
              <FontAwesomeIcon icon={faMinus} />
            </Button>
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
              <SvgContainer
                className="p-3"
                title={plotTitle(plots[index])}
                downloadName={plots[index].plotPath.split('/').slice(-1)[0]}
                plotPath={plots[index].plotURL}
                height="500px"
              />
            </div>
          )}
        </div>
      );
    }
    return display.slice(1);
  };

  function plotTitle(plot) {
    const { profileName, rsSet, signatureName } = plot;
    return `${profileName}/${rsSet}: Mutational Signature Profile of ${signatureName}`;
  }

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
        <Row className="msigportal-bg">
          <Col lg="auto">
            <CustomSelect
              id="mspSource"
              label="Signature Source"
              value={plots[0].signatureSource}
              options={plots[0].signatureSourceOptions}
              onChange={(source) => handleSource(source, 0)}
            />
          </Col>
          <Col lg="auto">
            <CustomSelect
              id="mspProfileName"
              label="Profile Name"
              value={plots[0].profileName}
              options={plots[0].profileNameOptions}
              onChange={(profile) => handleProfile(profile, 0)}
            />
          </Col>
          <Col lg="auto">
            <CustomSelect
              id="mspSet"
              label="Reference Signature Set"
              value={plots[0].rsSet}
              options={plots[0].rsSetOptions}
              onChange={(set) => handleSet(set, 0)}
            />
          </Col>
          <Col lg="auto">
            <CustomSelect
              id="mspStrategy"
              label="Experimental Strategy"
              value={plots[0].strategy}
              options={plots[0].strategyOptions}
              onChange={(strategy) => handleStrategy(strategy, 0)}
            />
          </Col>
          <Col lg="auto">
            <CustomSelect
              id="mspSigName"
              label="Signature Name"
              value={plots[0].signatureName}
              options={plots[0].signatureNameOptions}
              onChange={(name) => handleName(name, 0)}
            />
          </Col>
        </Row>
        <AdditionalControls />
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
      </Form>

      <div id="plot0">
        <div style={{ display: err ? 'block' : 'none' }}>
          <p>An error has occured. Please verify your input.</p>
        </div>
        {plots[0].plotURL && (
          <>
            <hr />
            <SvgContainer
              className="p-3"
              title={plotTitle(plots[0])}
              downloadName={plots[0].plotPath.split('/').slice(-1)[0]}
              plotPath={plots[0].plotURL}
              height="500px"
            />
          </>
        )}
      </div>
      {additionalPlots()}
    </div>
  );
}
