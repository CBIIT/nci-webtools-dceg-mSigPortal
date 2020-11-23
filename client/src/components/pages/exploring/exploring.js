import React, { useEffect } from 'react';
import { Nav } from 'react-bootstrap';
import { Route, NavLink, Redirect } from 'react-router-dom';
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
  const { publicDataOptions, displayTab } = useSelector(
    (state) => state.exploring
  );

  useEffect(() => {
    if (!Object.keys(publicDataOptions).length) populateControls();
  }, []);

  async function populateControls() {
    dispatchExpExposure({ loading: true });

    try {
      let [signatureData, exposureData] = await Promise.all([
        (await fetch(`api/public/Others/json/Exploring-Signature.json`)).json(),
        (await fetch('api/public/Others/json/Exploring-Exposure.json')).json(),
      ]);

      setInitalRefSigData(signatureData);
      const studyOptions = [...new Set(exposureData.map((data) => data.Study))];
      // default study
      const study = 'PCAWG';

      const cancerOptions = [
        ...new Set(
          exposureData
            .filter((data) => data.Study == study)
            .map((data) => data.Cancer_Type)
        ),
      ];
      // default cancer type
      const cancer = 'Lung-AdenoCA';

      const strategyOptions = [
        ...new Set(
          exposureData
            .filter((data) => data.Study == study)
            .map((data) => data.Dataset)
        ),
      ];
      const refSignatureSetOptions = [
        ...new Set(
          exposureData
            .filter((row) => row.Study == study)
            .map((row) => row.Signature_set_name)
        ),
      ];
      const refSignatureSet = refSignatureSetOptions[0];
      const signatureNameOptions = [
        ...new Set(
          signatureData
            .filter((row) => row.Signature_set_name == refSignatureSet)
            .map((row) => row.Signature_name)
        ),
      ];

      const params = {
        study: study,
        studyOptions: studyOptions,
        strategy: strategyOptions[0],
        strategyOptions: strategyOptions,
        cancer: cancer,
        cancerOptions: cancerOptions,
        refSigData: signatureData,
        refSignatureSet: refSignatureSet,
        refSignatureSetOptions: refSignatureSetOptions,
        signatureNameOptions: signatureNameOptions,
      };

      dispatchExploring({ publicDataOptions: exposureData });
      dispatchExpExposure(params);
    } catch (err) {
      dispatchError(err);
    }
    dispatchExpExposure({ loading: false });
  }

  async function setInitalRefSigData(data) {
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

  const links = [
    { name: 'Etiology', pathId: 'etiology' },
    { name: 'Signature', pathId: 'signature' },
    { name: 'Exposure', pathId: 'exposure' },
    { name: 'Download', pathId: 'download' },
  ];

  return (
    <div className="px-0">
      <div
        className="bg-white border border-top-0 mx-auto"
        style={{ width: 'max-content' }}
      >
        <Nav defaultActiveKey="summary">
          {links.map(({ name, pathId }) => (
            <div key={pathId} className="d-inline-block">
              <NavLink
                className="secondary-navlinks px-3 py-1 d-inline-block"
                activeClassName="active-secondary-navlinks"
                style={{
                  textDecoration: 'none',
                  fontSize: '11pt',
                  color: 'black',
                  fontWeight: '500',
                }}
                exact={true}
                to={`/exploring/${pathId}`}
              >
                {name}
              </NavLink>
            </div>
          ))}
        </Nav>
      </div>
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
            path="/exploring/exposure"
            render={() => (
              <ExposureExploring populateControls={populateControls} />
            )}
          />
          <Route path="/exploring/download" component={Download} />
        </div>
      </div>
    </div>
  );
}
