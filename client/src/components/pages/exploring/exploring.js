import React, { useEffect } from 'react';
import { Nav } from 'react-bootstrap';
import { Route, NavLink, Redirect } from 'react-router-dom';
import { useSelector } from 'react-redux';
import SignatureExploring from './signature/main';
import ExposureExploring from './exposure/main';
import EtiologyExploring from './etiology/main';
import {
  dispatchError,
  dispatchExploring,
  dispatchExpMutationalProfiles,
  dispatchExpCosineSimilarity,
  dispatchExpMutationalSigComparison,
  dispatchExpTumor,
  dispatchExpActivity,
  dispatchExpAssociation,
  dispatchExpDecomposition,
  dispatchExpLandscape,
  dispatchExpPrevalence,
  dispatchExpExposure,
} from '../../../services/store';
import './exploring.scss';

export default function Explore() {
  const rootURL = window.location.pathname;
  const {
    displayTab,
    projectID,
    refSigData,
    publicDataOptions,
    submitted,
  } = useSelector((state) => state.exploring);

  useEffect(() => {
    if (!Object.keys(refSigData).length) getInitalRefSigData();
    if (!Object.keys(publicDataOptions).length) getPublicDataOptions();
  }, []);

  async function getReferenceSignatureData(
    columns = [
      'Source',
      'Profile',
      'Signature_set_name',
      'Dataset',
      'Signature_name',
    ],
    filters = {}
  ) {
    try {
      const res = await fetch(`${rootURL}getReferenceSignatureData`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          args: {
            columns: columns,
            filters: filters,
          },
        }),
      });
      if (!res.ok) {
        // log error
      } else {
        return await res.json();
      }
    } catch (err) {}
  }

  async function getPublicDataOptions() {
    dispatchExpExposure({ loading: true });

    try {
      let [publicData, refSigData] = await Promise.all([
        fetch(`${rootURL}getPublicDataOptions`, {
          method: 'POST',
          headers: {
            Accept: 'application/json',
            'Content-Type': 'application/json',
          },
        }),
        getReferenceSignatureData(['Signature_set_name', 'Signature_name']),
      ]);

      if (!publicData.ok) {
        // console.log(await response.json());
      } else {
        const data = await publicData.json();
        const studyOptions = [...new Set(data.map((data) => data.Study))];
        // default study
        const study = 'PCAWG';

        const cancerOptions = [
          ...new Set(
            data
              .filter((data) => data.Study == study)
              .map((data) => data.Cancer_Type)
          ),
        ];
        // default cancer type
        const cancer = 'Lung-AdenoCA';

        const strategyOptions = [
          ...new Set(
            data
              .filter((data) => data.Study == study)
              .map((data) => data.Dataset)
          ),
        ];
        const refSignatureSetOptions = [
          ...new Set(
            refSigData.output.data.map((row) => row.Signature_set_name)
          ),
        ];
        const refSignatureSet = refSignatureSetOptions[0];
        const signatureNameOptions = [
          ...new Set(
            refSigData.output.data
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
          refSigData: refSigData.output.data,
          refSignatureSet: refSignatureSet,
          refSignatureSetOptions: refSignatureSetOptions,
          signatureNameOptions: signatureNameOptions,
        };

        dispatchExploring({ publicDataOptions: data });
        dispatchExpExposure(params);
      }
    } catch (err) {
      dispatchError(err);
    }
    dispatchExpExposure({ loading: false });
  }

  async function getInitalRefSigData() {
    // set loading indicators
    dispatchExpMutationalProfiles({ loading: true });
    dispatchExpCosineSimilarity({
      loading: true,
    });
    dispatchExpMutationalSigComparison({ loading: true });

    let data = (await getReferenceSignatureData()).output.data;
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
    { name: 'Etiology Exploring', pathId: 'etiology' },
    { name: 'Signature Exploring', pathId: 'signature' },
    { name: 'Exposure Exploring', pathId: 'exposure' },
  ];

  return (
    <div className="px-0">
      <div className="mx-3">
        <div className="mx-3 bg-white border border-top-0">
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
                <div className="d-md-none w-100"></div>
              </div>
            ))}
          </Nav>
        </div>
      </div>
      <div className="mx-3">
        <div className="mx-3 my-3 bg-white border">
          <div className="p-3">
            <Route
              exact
              path={`/exploring`}
              render={() => <Redirect to="/exploring/signature" />}
            />
            <Route path="/exploring/etiology" component={EtiologyExploring} />
            <Route path="/exploring/signature" component={SignatureExploring} />
            <Route path="/exploring/exposure" component={ExposureExploring} />
          </div>
        </div>
      </div>
    </div>
  );
}
