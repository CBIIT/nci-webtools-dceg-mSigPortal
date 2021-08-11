import React, { useEffect } from 'react';
import { Nav } from 'react-bootstrap';
import { Route, Redirect, NavLink } from 'react-router-dom';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../services/store/catalog';
import { actions as modalActions } from '../../../services/store/modal';
import { getJSON } from '../../../services/utils';
import Signature from './signature/index';
import Etiology from './etiology/etiology';

const actions = { ...catalogActions, ...modalActions };

export default function Explore() {
  const dispatch = useDispatch();
  const catalog = useSelector((state) => state.catalog);
  const { exposureSignature, displayTab } = catalog.catalog;
  const mergeCatalog = (state) =>
    dispatch(actions.mergeCatalog({ catalog: state }));
  const mergeSigMutationalProfiles = (state) =>
    dispatch(actions.mergeCatalog({ sigMutationalProfiles: state }));
  const mergeSigMutationalSigComparison = (state) =>
    dispatch(actions.mergeCatalog({ sigMutationalSigComparison: state }));
  const mergeSigCosineSimilarity = (state) =>
    dispatch(actions.mergeCatalog({ sigCosineSimilarity: state }));
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

    mergeCatalog({
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
                    className="secondary-navlinks px-3 py-1 d-inline-block catalogNav"
                    activeClassName="active-secondary-navlinks"
                    style={{
                      textDecoration: 'none',
                      fontSize: '12pt',
                      color: 'black',
                      fontWeight: '500',
                    }}
                    exact={true}
                    to={`/catalog/${id}`}
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
                    to={`/catalog/${id}`}
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
            path={`/catalog`}
            render={() => <Redirect to={`/catalog/${displayTab}`} />}
          />
          <Route path="/catalog/etiology" component={Etiology} />
          <Route path="/catalog/signature" component={Signature} />
        </div>
      </div>
    </div>
  );
}
