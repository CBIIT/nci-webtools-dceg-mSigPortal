import React, { useEffect } from 'react';
import { Tab, Nav } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
// import Overview from './overview';
// import Profile from './rsProfile';
import RsInMsigportal from './rsInMsigportal/rsInMsigportal-plot';
import RsProfile from './rsProfile/rsProfile';
import CosineSimilarity from './cosineSimilarity/cosineSimilarity';
import Comparison from './rsComparison';
import Download from './download';
import { getJSON } from '../../../../services/utils';
import { actions } from '../../../../services/store/catalog';

export default function ReferenceSignature() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.catalog);
  const mergeState = (state) =>
    dispatch(actions.mergeCatalog({ referenceSignature: state }));
  const mergeSigMutationalProfiles = (state) =>
    dispatch(actions.mergeCatalog({ sigMutationalProfiles: state }));
  const mergeSigMutationalSigComparison = (state) =>
    dispatch(actions.mergeCatalog({ sigMutationalSigComparison: state }));
  const mergeSigCosineSimilarity = (state) =>
    dispatch(actions.mergeCatalog({ sigCosineSimilarity: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const { projectID, displayTab, exposureSignature } = store.referenceSignature;

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

    mergeState({
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

  function submitR(fn, args) {
    return fetch(`web/explorationWrapper`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        fn: fn,
        args: args,
        projectID: projectID,
      }),
    });
  }

  const tabs = [
    // {
    //   component: <Overview submitR={(fn, args) => submitR(fn, args)} />,
    //   key: 'overview',
    //   title: 'RS in mSigPortal',
    // },
    {
      component: <RsInMsigportal />,
      key: 'overview',
      title: 'RS In mSigportal',
    },
    {
      component: <RsProfile />,
      key: 'RSProfile',
      title: 'RS Profile',
    },
    {
      component: <CosineSimilarity />,
      key: 'cosineSimilarity',
      title: 'RS Cosine Similarity',
    },
    {
      component: <Comparison submitR={(fn, args) => submitR(fn, args)} />,
      key: 'mutationalSignatureComparison',
      title: 'RS Comparison',
    },
    {
      component: <Download />,
      key: 'download',
      title: 'Download',
    },
  ];

  return (
    <div style={{ minHeight: '500px' }}>
      <Tab.Container
        transition={false}
        className="mt-2"
        activeKey={displayTab}
        onSelect={(tab) => mergeState({ displayTab: tab })}
      >
        <Nav variant="tabs">
          {tabs.map(({ key, title }) => (
            <Nav.Item key={key}>
              <Nav.Link
                eventKey={key}
                as="button"
                className="outline-none catalog"
              >
                <strong>{title}</strong>
              </Nav.Link>
            </Nav.Item>
          ))}
        </Nav>
        <Tab.Content
          className={`bg-white tab-pane-bordered rounded-0 d-block`}
          style={{ overflowX: 'auto' }}
        >
          {tabs.map(({ key, component }) => (
            <Tab.Pane key={key} eventKey={key} className="border-0">
              {component}
            </Tab.Pane>
          ))}
        </Tab.Content>
      </Tab.Container>
    </div>
  );
}
