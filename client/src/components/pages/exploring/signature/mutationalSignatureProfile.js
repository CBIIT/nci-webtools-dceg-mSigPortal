import React from 'react';
import { Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchExpMutationalProfiles,
} from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';
import Select from '../../../controls/select/select';

export default function MutationalSignatureProfile({ submitR }) {
  const rootURL = window.location.pathname;
  const { plots, debugR, err, loading } = useSelector(
    (state) => state.expMutationalProfiles
  );
  const { displayTab, refSigData } = useSelector((state) => state.exploring);

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
    dispatchExpMutationalProfiles({
      loading: false,
    });
    setRPlot(newPlots);

    // dispatchExpMutationalProfiles({
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
    //     dispatchExpMutationalProfiles({
    //       debugR: debugR + data.debugR,
    //     });
    //   });
    //   dispatchExpMutationalProfiles({
    //     loading: false,
    //   });
    //   setRPlot(newPlots);
    // } catch (err) {
    //   dispatchError(err);
    //   dispatchExpMutationalProfiles({ loading: false });
    // }
  }

  async function setRPlot(plots) {
    try {
      const svgs = await Promise.all(
        plots.map((plot) =>
          fetch(`${rootURL}public/${plot.plotPath}`).then((res) => res.blob())
        )
      );
      let newPlots = plots.slice();
      svgs.forEach((blob, index) => {
        if (newPlots[index].plotURL)
          URL.revokeObjectURL(newPlots[index].plotURL);
        const objectURL = URL.createObjectURL(blob);
        newPlots[index] = { ...newPlots[index], plotURL: objectURL };
      });

      dispatchExpMutationalProfiles({
        plots: newPlots,
      });
    } catch (err) {
      dispatchError(err);
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

    dispatchExpMutationalProfiles({
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

    dispatchExpMutationalProfiles({
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

    dispatchExpMutationalProfiles({
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

    dispatchExpMutationalProfiles({
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

    dispatchExpMutationalProfiles({
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

    dispatchExpMutationalProfiles({ plots: newPlots });
  }

  const additionalControls = () => {
    let controls = [];
    for (let index in plots) {
      controls.push(
        <Row className="justify-content-center" key={'control' + index}>
          <Col sm="2">
            <Select
              id={`mspSource${index}`}
              label="Signature Source"
              value={plots[index].signatureSource}
              options={plots[index].signatureSourceOptions}
              onChange={(source) => handleSource(source, index)}
            />
          </Col>
          <Col sm="2">
            <Select
              id={'mspProfileName' + index}
              label="Profile Name"
              value={plots[index].profileName}
              options={plots[index].profileNameOptions}
              onChange={(profile) => handleProfile(profile, index)}
            />
          </Col>
          <Col sm="3">
            <Select
              id={'mspSet' + index}
              label="Reference Signature Set"
              value={plots[index].refSignatureSet}
              options={plots[index].refSignatureSetOptions}
              onChange={(profile) => handleSet(profile, index)}
            />
          </Col>
          <Col sm="2">
            <Select
              id={'mspStrategy' + index}
              label="Experimental Strategy"
              value={plots[index].strategy}
              options={plots[index].strategyOptions}
              onChange={(strategy) => handleStrategy(strategy, index)}
            />
          </Col>
          <Col sm="2">
            <Select
              id={'mspSigName' + index}
              label="Signature Name"
              value={plots[index].signatureName}
              options={plots[index].signatureNameOptions}
              onChange={(name) => {
                let newPlots = plots.slice();
                newPlots[index] = { ...newPlots[index], signatureName: name };
                dispatchExpMutationalProfiles({
                  plots: newPlots,
                });
              }}
            />
          </Col>
          <Col sm="1" className="m-auto">
            <Button variant="secondary" onClick={() => removePlots(index)}>
              -
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
            <p>An error has occured. Check the debug section for more info.</p>
          </div>
          {plots[index].plotURL && (
            <div style={{ display: plots[index].plotURL ? 'block' : 'none' }}>
              <Plot
                plotName={plots[index].plotPath.split('/').slice(-1)[0]}
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
      <LoadingOverlay active={loading} />
      <Row className="justify-content-center">
        <Col sm="2">
          <Select
            id="mspSource"
            label="Signature Source"
            value={plots[0].signatureSource}
            options={plots[0].signatureSourceOptions}
            onChange={(source) => handleSource(source, 0)}
          />
        </Col>
        <Col sm="2">
          <Select
            id="mspProfileName"
            label="Profile Name"
            value={plots[0].profileName}
            options={plots[0].profileNameOptions}
            onChange={(profile) => handleProfile(profile, 0)}
          />
        </Col>
        <Col sm="3">
          <Select
            id="mspSet"
            label="Reference Signature Set"
            value={plots[0].refSignatureSet}
            options={plots[0].refSignatureSetOptions}
            onChange={(set) => handleSet(set, 0)}
          />
        </Col>
        <Col sm="2">
          <Select
            id="mspStrategy"
            label="Experimental Strategy"
            value={plots[0].strategy}
            options={plots[0].strategyOptions}
            onChange={(strategy) => handleStrategy(strategy, 0)}
          />
        </Col>
        <Col sm="2">
          <Select
            id="mspSigName"
            label="Signature Name"
            value={plots[0].signatureName}
            options={plots[0].signatureNameOptions}
            onChange={(name) => {
              let newPlots = plots.slice();
              newPlots[0] = { ...newPlots[0], signatureName: name };
              dispatchExpMutationalProfiles({
                plots: newPlots,
              });
            }}
          />
        </Col>
        <Col sm="1" className="m-auto">
          <Button variant="primary" onClick={handleCalculate}>
            Calculate
          </Button>
        </Col>
      </Row>
      {additionalControls()}
      <Row>
        <Col sm="1">
          <Button variant="secondary" onClick={() => addPlots()}>
            +
          </Button>
        </Col>
      </Row>
      <div id="plot0">
        <div style={{ display: err ? 'block' : 'none' }}>
          <p>An error has occured. Check the debug section for more info.</p>
        </div>
        <div style={{ display: plots[0].plotURL ? 'block' : 'none' }}>
          <Plot
            plotName={plots[0].plotPath.split('/').slice(-1)[0]}
            plotURL={plots[0].plotURL}
            maxHeight="1000px"
          />
        </div>
      </div>
      {additionalPlots()}
      <Debug msg={debugR} />
    </div>
  );
}
