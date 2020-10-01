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
  const {
    signatureSource,
    signatureSourceOptions,
    profileName,
    profileNameOptions,
    refSignatureSet,
    refSignatureSetOptions,
    strategy,
    strategyOptions,
    signatureName,
    signatureNameOptions,
    plotPath,
    plotURL,
    debugR,
    err,
    displayDebug,
    loading,
  } = useSelector((state) => state.expMutationalProfiles);
  const { displayTab, refSigData } = useSelector((state) => state.exploring);

  async function calculateR(fn, args) {
    dispatchExpMutationalProfiles({
      loading: true,
      err: false,
      debugR: '',
    });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();

        dispatchExpMutationalProfiles({
          loading: false,
          debugR: err,
        });
      } else {
        const { debugR, output } = await response.json();

        dispatchExpMutationalProfiles({
          debugR: debugR,
          loading: false,
          plotPath: output.plotPath,
        });
        setRPlot(output.plotPath, 'within');
      }
    } catch (err) {
      dispatchError(err);
      dispatchExpMutationalProfiles({ loading: false });
    }
  }

  async function setRPlot(plotPath) {
    if (plotPath) {
      try {
        const response = await fetch(`${rootURL}getPublicSVG`, {
          method: 'POST',
          headers: {
            Accept: 'image/svg',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({ path: plotPath }),
        });
        if (!response.ok) {
          // console.log(await response.json());
        } else {
          const pic = await response.blob();
          const objectURL = URL.createObjectURL(pic);

          if (plotURL) URL.revokeObjectURL(plotURL);
          dispatchExpMutationalProfiles({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchExpMutationalProfiles({ err: true, plotURL: '' });
    }
  }

  function handleSource(source) {
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

    dispatchExpMutationalProfiles({
      signatureSource: source,
      profileName: profileName,
      profileNameOptions: profileNameOptions,
      refSignatureSet: refSignatureSet,
      refSignatureSetOptions: refSignatureSetOptions,
      strategy: strategy,
      strategyOptions: strategyOptions,
      signatureName: signatureNameOptions[0],
      signatureNameOptions: signatureNameOptions,
    });
  }

  function handleProfile(profile) {
    const filteredData = refSigData.filter(
      (row) => row.Source == signatureSource && row.Profile == profile
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
              row.Source == signatureSource &&
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
              row.Source == signatureSource &&
              row.Profile == profile &&
              row.Signature_set_name == refSignatureSet &&
              row.Dataset == strategy
          )
          .map((row) => row.Signature_name)
      ),
    ];

    dispatchExpMutationalProfiles({
      profileName: profile,
      refSignatureSet: refSignatureSet,
      refSignatureSetOptions: refSignatureSetOptions,
      strategy: strategy,
      strategyOptions: strategyOptions,
      signatureName: signatureNameOptions[0],
      signatureNameOptions: signatureNameOptions,
    });
  }

  function handleSet(set) {
    let filteredData = refSigData.filter(
      (row) =>
        row.Source == signatureSource &&
        row.Profile == profileName &&
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
              row.Source == signatureSource &&
              row.Profile == profileName &&
              row.Signature_set_name == set &&
              row.Dataset == strategy
          )
          .map((row) => row.Signature_name)
      ),
    ];

    dispatchExpMutationalProfiles({
      refSignatureSet: set,
      strategy: strategy,
      strategyOptions: strategyOptions,
      signatureName: signatureNameOptions[0],
      signatureNameOptions: signatureNameOptions,
    });
  }

  function handleStrategy(strategy) {
    let filteredData = refSigData.filter(
      (row) =>
        row.Source == signatureSource &&
        row.Profile == profileName &&
        row.Signature_set_name == refSignatureSet &&
        row.Dataset == strategy
    );

    const signatureNameOptions = [
      ...new Set(filteredData.map((row) => row.Signature_name)),
    ];

    dispatchExpMutationalProfiles({
      strategy: strategy,
      signatureName: signatureNameOptions[0],
      signatureNameOptions: signatureNameOptions,
    });
  }

  return (
    <div>
      <LoadingOverlay active={loading} />
      <Row className="justify-content-center">
        <Col sm="2">
          <Select
            id="mspSource"
            label="Signature Source"
            value={signatureSource}
            options={signatureSourceOptions}
            onChange={handleSource}
          />
        </Col>
        <Col sm="2">
          <Select
            id="mspProfileName"
            label="Profile Name"
            value={profileName}
            options={profileNameOptions}
            onChange={handleProfile}
          />
        </Col>
        <Col sm="3">
          <Select
            id="mspSet"
            label="Reference Signature Set"
            value={refSignatureSet}
            options={refSignatureSetOptions}
            onChange={handleSet}
          />
        </Col>
        <Col sm="2">
          <Select
            id="mspStrategy"
            label="Experimental Strategy"
            value={strategy}
            options={strategyOptions}
            onChange={handleStrategy}
          />
        </Col>
        <Col sm="2">
          <Select
            id="mspSigName"
            label="Signature Name"
            value={signatureName}
            options={signatureNameOptions}
            onChange={(name) =>
              dispatchExpMutationalProfiles({
                signatureName: name,
              })
            }
          />
        </Col>
        <Col sm="1" className="m-auto">
          <Button
            variant="primary"
            onClick={() => {
              calculateR('mutationalProfiles', {
                signatureSource: signatureSource,
                profileName: profileName,
                refSignatureSet: refSignatureSet,
                experimentalStrategy: strategy,
                signatureName: signatureName,
              });
            }}
          >
            Calculate
          </Button>
        </Col>
      </Row>

      <div id="withinPlot">
        <div style={{ display: err ? 'block' : 'none' }}>
          <p>An error has occured. Check the debug section for more info.</p>
        </div>
        <div style={{ display: plotURL ? 'block' : 'none' }}>
          <Plot
            plotName={plotPath.split('/').slice(-1)[0]}
            plotURL={plotURL}
            maxHeight="1000px"
          />
        </div>
      </div>

      <Debug msg={debugR} />
    </div>
  );
}
