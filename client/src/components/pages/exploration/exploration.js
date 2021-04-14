import React, { useEffect } from 'react';
import { Nav } from 'react-bootstrap';
import { Route, Redirect, NavLink } from 'react-router-dom';
import { useSelector, useDispatch } from 'react-redux';
import {
  actions as explorationActions,
  getInitialState,
} from '../../../services/store/exploration';
import { actions as modalActions } from '../../../services/store/modal';
import Signature from './signature/index';
import Exposure from './exposure';
import Aetiology from './aetiology';
import Download from './download';
import './exploration.scss';

const actions = { ...explorationActions, ...modalActions };

export default function Explore() {
  const dispatch = useDispatch();
  const exploration = useSelector((state) => state.exploration);
  const { exposureSignature, displayTab } = exploration.exploration;
  const mergeExploration = (state) =>
    dispatch(actions.mergeExploration({ exploration: state }));
  const mergeExposure = (state) =>
    dispatch(actions.mergeExploration({ exposure: state }));
  const mergeSigMutationalProfiles = (state) =>
    dispatch(actions.mergeExploration({ sigMutationalProfiles: state }));
  const mergeSigMutationalSigComparison = (state) =>
    dispatch(actions.mergeExploration({ sigMutationalSigComparison: state }));
  const mergeSigCosineSimilarity = (state) =>
    dispatch(actions.mergeExploration({ sigCosineSimilarity: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  useEffect(() => {
    if (!exposureSignature.length) populateControls();
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
      let [
        signatureData,
        exposureCancer,
        exposureSignature,
        signatureNames,
      ] = await Promise.all([
        getFileS3(`msigportal/Database/Others/json/Exploring-Signature.json`),
        getFileS3(
          'msigportal/Database/Others/json/Exploring-Exposure-cancertype.json'
        ),
        getFileS3('msigportal/Database/Others/json/Exploring-Exposure.json'),
        getFileS3('msigportal/Database/Others/json/Signature_name.json'),
      ]);

      populateSignatureExp(signatureData);
      populateExposureExp(exposureCancer, exposureSignature, signatureNames);
    } catch (err) {
      mergeError(err.message);
    }
  }

  async function populateExposureExp(
    exposureCancer,
    exposureSignature,
    signatureNames
  ) {
    mergeExposure({ loading: true });

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

    mergeExposure(params);

    mergeExploration({
      exposureCancer: exposureCancer,
      exposureSignature: exposureSignature,
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

    mergeExploration({
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
                    className="secondary-navlinks px-3 py-1 d-inline-block explorationNav"
                    activeClassName="active-secondary-navlinks"
                    style={{
                      textDecoration: 'none',
                      fontSize: '11pt',
                      color: 'black',
                      fontWeight: '500',
                    }}
                    exact={true}
                    to={`/exploration/${id}`}
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
                    to={`/exploration/${id}`}
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
            path={`/exploration`}
            render={() => <Redirect to={`/exploration/${displayTab}`} />}
          />
          <Route path="/exploration/aetiology" component={Aetiology} />
          <Route path="/exploration/signature" component={Signature} />
          <Route
            path="/exploration/exposure/:exampleName?"
            render={(props) => (
              <Exposure {...props} populateControls={populateControls} />
            )}
          />
          <Route path="/exploration/download" component={Download} />
        </div>
      </div>
    </div>
  );
}
