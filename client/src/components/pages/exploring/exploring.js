import React, { useEffect } from 'react';
import { Route, Redirect } from 'react-router-dom';
import { useSelector } from 'react-redux';
import SignatureExploring from './signature';
import ExposureExploring from './exposure';
import EtiologyExploring from './etiology';
import Download from './download';
import {
  dispatchError,
  dispatchExploring,
  dispatchExpMutationalProfiles,
  dispatchExpCosineSimilarity,
  dispatchExpMutationalSigComparison,
  dispatchExpExposure,
} from '../../../services/store';
import './exploring.scss';

export default function Explore() {
  const { exposureSignature, displayTab } = useSelector(
    (state) => state.exploring
  );

  useEffect(() => {
    if (!exposureSignature.length) populateControls();
  }, []);

  async function populateControls() {
    try {
      let [
        signatureData,
        exposureCancer,
        exposureSignature,
        signatureNames,
      ] = await Promise.all([
        (await fetch(`api/public/Others/json/Exploring-Signature.json`)).json(),
        (
          await fetch(
            'api/public/Others/json/Exploring-Exposure-cancertype.json'
          )
        ).json(),
        (
          await fetch(
            'api/public/Others/json/Exploring-Exposure-Signature-set-name.json'
          )
        ).json(),
        (await fetch('api/public/Others/json/Signature_name.json')).json(),
      ]);

      populateSignatureExp(signatureData);
      populateExposureExp(exposureCancer, exposureSignature, signatureNames);
    } catch (err) {
      dispatchError(err);
    }
  }

  async function populateExposureExp(
    exposureCancer,
    exposureSignature,
    signatureNames
  ) {
    dispatchExpExposure({ loading: true });

    const studyOptions = [
      ...new Set(exposureSignature.map((data) => data.Study)),
    ];
    const study = 'PCAWG'; // default

    const strategyOptions = [
      ...new Set(
        exposureSignature
          .filter((data) => data.Study == study)
          .map((data) => data.Dataset)
      ),
    ];
    const strategy = strategyOptions[0];

    const refSignatureSetOptions = [
      ...new Set(
        exposureSignature
          .filter((row) => row.Study == study && row.Dataset == strategy)
          .map((row) => row.Signature_set_name)
      ),
    ];
    const refSignatureSet = 'COSMIC v3 Signatures (SBS)'; // default

    const cancerOptions = [
      ...new Set(
        exposureCancer
          .filter((data) => data.Study == study && data.Dataset == strategy)
          .map((data) => data.Cancer_Type)
      ),
    ];
    const cancer = 'Lung-AdenoCA'; // default

    const signatureNameOptions = [
      ...new Set(
        signatureNames
          .filter((row) => row.Signature_set_name == refSignatureSet)
          .map((row) => row.Signature_name)
      ),
    ];

    const params = {
      study: study,
      studyOptions: studyOptions,
      strategy: strategy,
      strategyOptions: strategyOptions,
      cancer: cancer,
      cancerOptions: cancerOptions,
      signatureNames: signatureNames,
      refSignatureSet: refSignatureSet,
      refSignatureSetOptions: refSignatureSetOptions,
      signatureNameOptions: signatureNameOptions,
      loading: false,
    };

    dispatchExpExposure(params);

    dispatchExploring({
      exposureCancer: exposureCancer,
      exposureSignature: exposureSignature,
      signatureNames: signatureNames,
    });
  }

  async function populateSignatureExp(data) {
    // set loading indicators
    dispatchExpMutationalProfiles({ loading: true });
    dispatchExpCosineSimilarity({
      loading: true,
    });
    dispatchExpMutationalSigComparison({ loading: true });

    const signatureSourceOptions = [...new Set(data.map((row) => row.Source))];
    const signatureSource = signatureSourceOptions[0];
    const profileNameOptions = [
      ...new Set(
        data
          .filter((row) => row.Source == signatureSource)
          .map((row) => row.Profile)
      ),
    ];
    const profileName = profileNameOptions[0];
    const refSignatureSetOptions = [
      ...new Set(
        data
          .filter(
            (row) => row.Source == signatureSource && row.Profile == profileName
          )
          .map((row) => row.Signature_set_name)
      ),
    ];
    const refSignatureSet = refSignatureSetOptions[0];
    const refSignatureSet2 =
      refSignatureSetOptions[1] || refSignatureSetOptions[0];

    const strategyOptions = [
      ...new Set(
        data
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
        data
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
    const signatureNameOptions2 = [
      ...new Set(
        data
          .filter(
            (row) =>
              row.Source == signatureSource &&
              row.Profile == profileName &&
              row.Signature_set_name == refSignatureSet2 &&
              row.Dataset == strategy
          )
          .map((row) => row.Signature_name)
      ),
    ];

    dispatchExploring({
      refSigData: data,
    });

    dispatchExpMutationalProfiles({
      plots: [
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
      loading: false,
    });

    dispatchExpCosineSimilarity({
      profileName: profileName,
      profileNameOptions: profileNameOptions,
      refSignatureSet1: refSignatureSetOptions[0],
      refSignatureSet2: refSignatureSetOptions[1] || refSignatureSetOptions[0],
      refSignatureSetOptions1: refSignatureSetOptions,
      refSignatureSetOptions2: refSignatureSetOptions,
      loading: false,
    });

    dispatchExpMutationalSigComparison({
      profileName: profileName,
      profileNameOptions: profileNameOptions,
      refSignatureSet1: refSignatureSet,
      refSignatureSet2: refSignatureSet2,
      refSignatureSetOptions1: refSignatureSetOptions,
      refSignatureSetOptions2: refSignatureSetOptions,
      signatureName1: signatureNameOptions[0],
      signatureNameOptions1: signatureNameOptions,
      signatureName2: signatureNameOptions2[0],
      signatureNameOptions2: signatureNameOptions2,
      loading: false,
    });
  }

  // const links = [
  //   { name: 'Etiology', pathId: 'etiology' },
  //   { name: 'Signature', pathId: 'signature' },
  //   { name: 'Exposure', pathId: 'exposure' },
  //   { name: 'Download', pathId: 'download' },
  // ];

  return (
    <div className="px-0">
      <div className="mx-3">
        <div className="mx-3 my-3">
          <Route
            exact
            path={`/exploring`}
            render={() => <Redirect to={`/exploring/${displayTab}`} />}
          />
          <Route path="/exploring/etiology" component={EtiologyExploring} />
          <Route path="/exploring/signature" component={SignatureExploring} />
          <Route
            path="/exploring/exposure/:exampleName?"
            render={(props) => (
              <ExposureExploring
                {...props}
                populateControls={populateControls}
              />
            )}
          />
          <Route path="/exploring/download" component={Download} />
        </div>
      </div>
    </div>
  );
}
