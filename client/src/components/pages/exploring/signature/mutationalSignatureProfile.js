import React from 'react';
import { Row, Col, Button, Form } from 'react-bootstrap';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';
import Select from '../../../controls/select/select';
import { useSelector, useDispatch } from 'react-redux';
import { actions as exploringActions } from '../../../../services/store/exploring';
import { actions as modalActions } from '../../../../services/store/modal';

const actions = { ...exploringActions, ...modalActions };

export default function MutationalSignatureProfile({ submitR }) {
  const dispatch = useDispatch();
  const exploring = useSelector((state) => state.exploring);
  const mergeSigMutationalProfiles = (state) =>
    dispatch(actions.mergeExploring({ sigMutationalProfiles: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const { plots, debugR, err, loading } = exploring.sigMutationalProfiles;
  const { displayTab, refSigData } = exploring.exploring;

  // dataFolder/Reference_Signature_Profiles_SVG/refSignatureSet/profileName+signatureName
  function buildPlotPath(profileName, refSignatureSet, signatureName) {
    const profile = profileName.match(/[a-z]+|\d+/gi).join('_');
    const set = refSignatureSet
      .replace(/\s/g, '_')
      .replace(/[^a-zA-Z0-9-_]/gi, '');

    return `Signature/Reference_Signature_Profiles_SVG/${set}/${profile}_plots_mSigPortal_${signatureName}.svg`;
  }

  async function handleCalculate() {
    const plotPaths = plots.map((plot, index) =>
      buildPlotPath(plot.profileName, plot.refSignatureSet, plot.signatureName)
    );

    let newPlots = plots.slice();
    plotPaths.forEach((path, index) => {
      newPlots[index] = {
        ...newPlots[index],
        plotPath: path,
      };
    });
    mergeSigMutationalProfiles({
      loading: false,
    });
    setRPlot(newPlots);

    // mergeSigMutationalProfiles({
    //   loading: true,
    //   err: false,
    //   debugR: '',
    // });
    // try {
    //   const reqR = await Promise.all(
    //     plots.map((plot, index) =>
    //       submitR('mutationalProfiles', {
    //         signatureSource: plot.signatureSource,
    //         profileName: plot.profileName,
    //         refSignatureSet: plot.refSignatureSet,
    //         experimentalStrategy: plot.strategy,
    //         signatureName: plot.signatureName,
    //       }).then((res) => res.json())
    //     )
    //   );
    //   let newPlots = plots.slice();
    //   reqR.forEach((data, index) => {
    //     newPlots[index] = {
    //       ...newPlots[index],
    //       plotPath: data.output.plotPath,
    //     };
    //     mergeSigMutationalProfiles({
    //       debugR: debugR + data.debugR,
    //     });
    //   });
    //   mergeSigMutationalProfiles({
    //     loading: false,
    //   });
    //   setRPlot(newPlots);
    // } catch (err) {
    //   mergeError(err.message );
    //   mergeSigMutationalProfiles({ loading: false });
    // }
  }

  async function setRPlot(plots) {
    try {
      const svgs = await Promise.all(
        plots.map((plot) =>
          fetch(`api/public/${plot.plotPath}`).then((res) => res.blob())
        )
      );
      let newPlots = plots.slice();
      svgs.forEach((blob, index) => {
        if (newPlots[index].plotURL)
          URL.revokeObjectURL(newPlots[index].plotURL);
        const objectURL = URL.createObjectURL(blob);
        newPlots[index] = { ...newPlots[index], plotURL: objectURL };
      });

      mergeSigMutationalProfiles({
        plots: newPlots,
      });
    } catch (err) {
      mergeError(err.message);
    }
  }

  function handleSource(source, index) {
    const filteredData = refSigData.filter((row) => row.Source == source);
    const profileNameOptions = [
      ...new Set(filteredData.map((row) => row.Profile)),
    ];
    const profileName = profileNameOptions[0];
    const refSignatureSetOptions = [
      ...new Set(
        filteredData
          .filter((row) => row.Source == source && row.Profile == profileName)
          .map((row) => row.Signature_set_name)
      ),
    ];
    const refSignatureSet = refSignatureSetOptions[0];
    const strategyOptions = [
      ...new Set(
        filteredData
          .filter(
            (row) =>
              row.Source == source &&
              row.Profile == profileName &&
              row.Signature_set_name == refSignatureSet
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
              row.Signature_set_name == refSignatureSet &&
              row.Dataset == strategy
          )
          .map((row) => row.Signature_name)
      ),
    ];

    let newPlots = plots.slice();
    newPlots[index] = {
      ...newPlots[index],
      signatureSource: source,
      profileName: profileName,
      profileNameOptions: profileNameOptions,
      refSignatureSet: refSignatureSet,
      refSignatureSetOptions: refSignatureSetOptions,
      strategy: strategy,
      strategyOptions: strategyOptions,
      signatureName: signatureNameOptions[0],
      signatureNameOptions: signatureNameOptions,
    };

    mergeSigMutationalProfiles({
      plots: newPlots,
    });
  }

  function handleProfile(profile, index) {
    const filteredData = refSigData.filter(
      (row) =>
        row.Source == plots[index].signatureSource && row.Profile == profile
    );
    const refSignatureSetOptions = [
      ...new Set(filteredData.map((row) => row.Signature_set_name)),
    ];
    const refSignatureSet = refSignatureSetOptions[0];
    const strategyOptions = [
      ...new Set(
        filteredData
          .filter(
            (row) =>
              row.Source == plots[index].signatureSource &&
              row.Profile == profile &&
              row.Signature_set_name == refSignatureSet
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
              row.Signature_set_name == refSignatureSet &&
              row.Dataset == strategy
          )
          .map((row) => row.Signature_name)
      ),
    ];

    let newPlots = plots.slice();
    newPlots[index] = {
      ...newPlots[index],
      profileName: profile,
      refSignatureSet: refSignatureSet,
      refSignatureSetOptions: refSignatureSetOptions,
      strategy: strategy,
      strategyOptions: strategyOptions,
      signatureName: signatureNameOptions[0],
      signatureNameOptions: signatureNameOptions,
    };

    mergeSigMutationalProfiles({
      plots: newPlots,
    });
  }

  function handleSet(set, index) {
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

    let newPlots = plots.slice();
    newPlots[index] = {
      ...newPlots[index],
      refSignatureSet: set,
      strategy: strategy,
      strategyOptions: strategyOptions,
      signatureName: signatureNameOptions[0],
      signatureNameOptions: signatureNameOptions,
    };

    mergeSigMutationalProfiles({
      plots: newPlots,
    });
  }

  function handleStrategy(strategy, index) {
    let filteredData = refSigData.filter(
      (row) =>
        row.Source == plots[index].signatureSource &&
        row.Profile == plots[index].profileName &&
        row.Signature_set_name == plots[index].refSignatureSet &&
        row.Dataset == strategy
    );

    const signatureNameOptions = [
      ...new Set(filteredData.map((row) => row.Signature_name)),
    ];

    let newPlots = plots.slice();
    newPlots[index] = {
      ...newPlots[index],
      strategy: strategy,
      signatureName: signatureNameOptions[0],
      signatureNameOptions: signatureNameOptions,
    };

    mergeSigMutationalProfiles({
      plots: newPlots,
    });
  }

  function addPlots() {
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
    const refSignatureSetOptions = [
      ...new Set(
        refSigData
          .filter(
            (row) => row.Source == signatureSource && row.Profile == profileName
          )
          .map((row) => row.Signature_set_name)
      ),
    ];
    const refSignatureSet = refSignatureSetOptions[0];

    const strategyOptions = [
      ...new Set(
        refSigData
          .filter(
            (row) =>
              row.Source == signatureSource &&
              row.Profile == profileName &&
              row.Signature_set_name == refSignatureSet
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
              row.Signature_set_name == refSignatureSet &&
              row.Dataset == strategy
          )
          .map((row) => row.Signature_name)
      ),
    ];

    mergeSigMutationalProfiles({
      plots: [
        ...plots,
        {
          signatureSource: signatureSource,
          signatureSourceOptions: signatureSourceOptions,
          profileName: profileName,
          profileNameOptions: profileNameOptions,
          refSignatureSet: refSignatureSet,
          refSignatureSetOptions: refSignatureSetOptions,
          strategy: strategy,
          strategyOptions: strategyOptions,
          signatureName: signatureNameOptions[0],
          signatureNameOptions: signatureNameOptions,
          plotPath: '',
          plotURL: '',
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
          <Col lg="2">
            <Select
              id={`mspSource${index}`}
              label="Signature Source"
              value={plots[index].signatureSource}
              options={plots[index].signatureSourceOptions}
              onChange={(source) => handleSource(source, index)}
            />
          </Col>
          <Col lg="2">
            <Select
              id={'mspProfileName' + index}
              label="Profile Name"
              value={plots[index].profileName}
              options={plots[index].profileNameOptions}
              onChange={(profile) => handleProfile(profile, index)}
            />
          </Col>
          <Col lg="3">
            <Select
              id={'mspSet' + index}
              label="Reference Signature Set"
              value={plots[index].refSignatureSet}
              options={plots[index].refSignatureSetOptions}
              onChange={(profile) => handleSet(profile, index)}
            />
          </Col>
          <Col lg="2">
            <Select
              id={'mspStrategy' + index}
              label="Experimental Strategy"
              value={plots[index].strategy}
              options={plots[index].strategyOptions}
              onChange={(strategy) => handleStrategy(strategy, index)}
            />
          </Col>
          <Col lg="2">
            <Select
              id={'mspSigName' + index}
              label="Signature Name"
              value={plots[index].signatureName}
              options={plots[index].signatureNameOptions}
              onChange={(name) => {
                let newPlots = plots.slice();
                newPlots[index] = { ...newPlots[index], signatureName: name };
                mergeSigMutationalProfiles({
                  plots: newPlots,
                });
              }}
            />
          </Col>
          <Col lg="1" className="d-flex">
            <Button
              className="my-auto text-danger"
              variant="link"
              onClick={() => removePlots(index)}
              title={'Remove Plot ' + (parseInt(index) + 1)}
              style={{ width: '2.25rem' }}
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
              <Plot
                className="p-3"
                title={plots[index].plotPath.split('/').slice(-1)[0]}
                downloadName={plots[index].plotPath.split('/').slice(-1)[0]}
                plotURL={plots[index].plotURL}
                maxHeight="1000px"
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
      <Form className="p-3">
        <LoadingOverlay active={loading} />
        <Row>
          <Col lg="2">
            <Select
              id="mspSource"
              label="Signature Source"
              value={plots[0].signatureSource}
              options={plots[0].signatureSourceOptions}
              onChange={(source) => handleSource(source, 0)}
            />
          </Col>
          <Col lg="2">
            <Select
              id="mspProfileName"
              label="Profile Name"
              value={plots[0].profileName}
              options={plots[0].profileNameOptions}
              onChange={(profile) => handleProfile(profile, 0)}
            />
          </Col>
          <Col lg="3">
            <Select
              id="mspSet"
              label="Reference Signature Set"
              value={plots[0].refSignatureSet}
              options={plots[0].refSignatureSetOptions}
              onChange={(set) => handleSet(set, 0)}
            />
          </Col>
          <Col lg="2">
            <Select
              id="mspStrategy"
              label="Experimental Strategy"
              value={plots[0].strategy}
              options={plots[0].strategyOptions}
              onChange={(strategy) => handleStrategy(strategy, 0)}
            />
          </Col>
          <Col lg="2">
            <Select
              id="mspSigName"
              label="Signature Name"
              value={plots[0].signatureName}
              options={plots[0].signatureNameOptions}
              onChange={(name) => {
                let newPlots = plots.slice();
                newPlots[0] = { ...newPlots[0], signatureName: name };
                mergeSigMutationalProfiles({
                  plots: newPlots,
                });
              }}
            />
          </Col>
          <Col />
        </Row>
        <AdditionalControls />
        <Row className="mt-3">
          <Col />
          <Col lg="1" className="d-flex">
            <Button
              className="my-auto"
              variant="link"
              onClick={() => addPlots()}
              title="Add Plot"
              style={{ width: '2.25rem', textDecoration: 'none' }}
            >
              <span className="text-nowrap" title="Add Plot">
                <FontAwesomeIcon icon={faPlus} /> Add Plot
              </span>
            </Button>
          </Col>
          <Col lg="1" className="d-flex">
            <Button
              className="my-auto"
              variant="primary"
              onClick={handleCalculate}
            >
              Calculate
            </Button>
          </Col>
        </Row>
      </Form>

      <div id="plot0">
        <div style={{ display: err ? 'block' : 'none' }}>
          <p>An error has occured. Please verify your input.</p>
        </div>
        <div style={{ display: plots[0].plotURL ? 'block' : 'none' }}>
          <hr />
          <Plot
            className="p-3"
            title={plots[0].plotPath.split('/').slice(-1)[0]}
            downloadName={plots[0].plotPath.split('/').slice(-1)[0]}
            plotURL={plots[0].plotURL}
            maxHeight="1000px"
          />
        </div>
      </div>
      {additionalPlots()}
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
