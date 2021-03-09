import React, { useEffect } from 'react';
import { Nav } from 'react-bootstrap';
import { Route, Redirect, NavLink } from 'react-router-dom';
import { useSelector, useDispatch } from 'react-redux';
import {
  actions as exploringActions,
  getInitialState,
} from '../../../services/store/exploring';
import { actions as modalActions } from '../../../services/store/modal';
import Signature from './signature/index';
import Exposure from './exposure';
import Aetiology from './aetiology';
import Download from './download';
import './exploring.scss';

const actions = { ...exploringActions, ...modalActions };

export default function Explore() {
  const dispatch = useDispatch();
  const exploring = useSelector((state) => state.exploring);
  const { exposureData, displayTab } = exploring.exploring;
  const mergeExploring = (state) =>
    dispatch(actions.mergeExploring({ exploring: state }));
  const mergeExposure = (state) =>
    dispatch(actions.mergeExploring({ exposure: state }));
  const mergeSigMutationalProfiles = (state) =>
    dispatch(actions.mergeExploring({ sigMutationalProfiles: state }));
  const mergeSigMutationalSigComparison = (state) =>
    dispatch(actions.mergeExploring({ sigMutationalSigComparison: state }));
  const mergeSigCosineSimilarity = (state) =>
    dispatch(actions.mergeExploring({ sigCosineSimilarity: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  useEffect(() => {
    if (!exposureData.length) populateControls();
  }, []);

  const getFileS3 = (path) =>
    fetch(`api/getFileS3`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        path: path,
      }),
    })
      .then((res) => res.json())
      .then((data) => data);

  async function populateControls() {
    try {
      let [signatureData, exposureData, signatureNames] = await Promise.all([
        getFileS3(`msigportal/Database/Others/json/Exploring-Signature.json`),
        getFileS3('msigportal/Database/Others/json/Exploring-Exposure.json'),

        getFileS3('msigportal/Database/Others/json/Signature_name.json'),
      ]);

      populateSignatureExp(signatureData);
      populateExposureExp(exposureData, signatureNames);
    } catch (err) {
      mergeError(err.message);
    }
  }

  async function populateExposureExp(exposureData, signatureNames) {
    mergeExposure({ loading: true });

    const studyOptions = [...new Set(exposureData.map((data) => data.Study))];
    const study = 'PCAWG'; // default

    const strategyOptions = [
      ...new Set(
        exposureData
          .filter((data) => data.Study == study)
          .map((data) => data.Dataset)
      ),
    ];
    const strategy = strategyOptions[0];

    const refSignatureSetOptions = [
      ...new Set(
        exposureData
          .filter((row) => row.Study == study && row.Dataset == strategy)
          .map((row) => row.Signature_set_name)
      ),
    ];
    const refSignatureSet = 'COSMIC v3 Signatures (SBS)'; // default

    const cancerOptions = [
      ...new Set(
        exposureData
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

    mergeExposure(params);

    mergeExploring({
      exposureData: exposureData,
      signatureNames: signatureNames,
    });
  }

  async function populateSignatureExp(data) {
    // set loading indicators
    mergeSigMutationalProfiles({ loading: true });
    mergeSigCosineSimilarity({
      loading: true,
    });
    mergeSigMutationalSigComparison({ loading: true });

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

    mergeExploring({
      refSigData: data,
    });

    mergeSigMutationalProfiles({
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

    mergeSigCosineSimilarity({
      profileName: profileName,
      profileNameOptions: profileNameOptions,
      refSignatureSet1: refSignatureSetOptions[0],
      refSignatureSet2: refSignatureSetOptions[1] || refSignatureSetOptions[0],
      refSignatureSetOptions1: refSignatureSetOptions,
      refSignatureSetOptions2: refSignatureSetOptions,
      loading: false,
    });

    mergeSigMutationalSigComparison({
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

  const tabs = [
    { name: 'Aetiology', id: 'aetiology' },
    { name: 'Signature', id: 'signature' },
    { name: 'Exposure', id: 'exposure' },
    { name: 'Download', id: 'download' },
  ];

  return (
    <div className="px-0">
      <div className="mx-3">
        <div className="mx-3 bg-white border border-top-0">
          {/* for desktops and tablets */}
          <div className="d-none d-md-block">
            <Nav defaultActiveKey="profilerSummary">
              {tabs.map(({ name, id }) => (
                <div key={id} className="d-inline-block">
                  <NavLink
                    className="secondary-navlinks px-3 py-1 d-inline-block exploringNav"
                    activeClassName="active-secondary-navlinks"
                    style={{
                      textDecoration: 'none',
                      fontSize: '11pt',
                      color: 'black',
                      fontWeight: '500',
                    }}
                    exact={true}
                    to={`/exploring/${id}`}
                  >
                    {name}
                  </NavLink>
                  <div className="d-md-none w-100"></div>
                </div>
              ))}
            </Nav>
          </div>
          {/* for mobile devices */}
          <div className="row d-md-none">
            <Nav defaultActiveKey="summary">
              {tabs.map(({ name, id }) => (
                <div key={id} className="col-12 text-center">
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
                    to={`/exploring/${id}`}
                  >
                    {name}
                  </NavLink>
                  <div className="d-md-none w-100"></div>
                </div>
              ))}
            </Nav>
          </div>
        </div>
        <div className="mx-3 my-3">
          <Route
            exact
            path={`/exploring`}
            render={() => <Redirect to={`/exploring/${displayTab}`} />}
          />
          <Route path="/exploring/aetiology" component={Aetiology} />
          <Route path="/exploring/signature" component={Signature} />
          <Route
            path="/exploring/exposure/:exampleName?"
            render={(props) => (
              <Exposure {...props} populateControls={populateControls} />
            )}
          />
          <Route path="/exploring/download" component={Download} />
        </div>
      </div>
    </div>
  );
}
