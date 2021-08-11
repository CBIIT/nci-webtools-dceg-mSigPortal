import React, { useEffect } from 'react';
import { Nav } from 'react-bootstrap';
import { Route, Redirect, NavLink } from 'react-router-dom';
import { useSelector, useDispatch } from 'react-redux';
import { actions as explorationActions } from '../../../services/store/exploration';
import { actions as modalActions } from '../../../services/store/modal';
import { getJSON } from '../../../services/utils';
import Signature from './signature/index';
import Etiology from './etiology/etiology';
import './exploration.scss';

const actions = { ...explorationActions, ...modalActions };

export default function Explore() {
  const dispatch = useDispatch();
  const exploration = useSelector((state) => state.exploration);
  const { exposureSignature, displayTab } = exploration.exploration;
  const mergeExploration = (state) =>
    dispatch(actions.mergeExploration({ exploration: state }));
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

  async function populateControls() {
    try {
      let signatureData = await getJSON(`Others/json/Exploring-Signature.json`);
      populateSignatureExp(signatureData);
    } catch (err) {
      mergeError(err.message);
    }
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
    const rsSetOptions = [
      ...new Set(
        data
          .filter(
            (row) => row.Source == signatureSource && row.Profile == profileName
          )
          .map((row) => row.Signature_set_name)
      ),
    ];
    const rsSet = rsSetOptions[0];
    const rsSet2 = rsSetOptions[1] || rsSetOptions[0];

    const strategyOptions = [
      ...new Set(
        data
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
        data
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
    const signatureNameOptions2 = [
      ...new Set(
        data
          .filter(
            (row) =>
              row.Source == signatureSource &&
              row.Profile == profileName &&
              row.Signature_set_name == rsSet2 &&
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
          rsSet: rsSet,
          rsSetOptions: rsSetOptions,
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
      rsSet1: rsSetOptions[0],
      rsSet2: rsSetOptions[1] || rsSetOptions[0],
      rsSetOptions1: rsSetOptions,
      rsSetOptions2: rsSetOptions,
      loading: false,
    });

    mergeSigMutationalSigComparison({
      profileName: profileName,
      profileNameOptions: profileNameOptions,
      rsSet1: rsSet,
      rsSet2: rsSet2,
      rsSetOptions1: rsSetOptions,
      rsSetOptions2: rsSetOptions,
      signatureName1: signatureNameOptions[0],
      signatureNameOptions1: signatureNameOptions,
      signatureName2: signatureNameOptions2[0],
      signatureNameOptions2: signatureNameOptions2,
      loading: false,
    });
  }

  const tabs = [
    { name: 'Etiology', id: 'etiology' },
    { name: 'Signature', id: 'signature' },
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
                      fontSize: '12pt',
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
                      fontSize: '12pt',
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
          <Route path="/exploration/etiology" component={Etiology} />
          <Route path="/exploration/signature" component={Signature} />
        </div>
      </div>
    </div>
  );
}
