import React, { useEffect } from 'react';
import { Card, Nav, Form } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import SignatureExploring from './signature/main';
import ExposureExploring from './exposure/main';
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

const { Header, Body } = Card;
const { Item, Link } = Nav;
const { Group, Label, Check } = Form;

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
          refSigData: refSigData.output.data,
          refSignatureSet: refSignatureSet,
          refSignatureSetOptions: refSignatureSetOptions,
          signatureNameOptions: signatureNameOptions,
        };
        const associationParams = {
          cancer: 'None',
          cancerOptions: ['None', ...cancerOptions],
        };
        const landscapeParams = {
          cancer: cancer,
          cancerOptions: cancerOptions,
        };

        dispatchExploring({ publicDataOptions: data });
        dispatchExpExposure(params);
        dispatchExpAssociation(associationParams);
        dispatchExpLandscape(landscapeParams);
        dispatchExpPrevalence(landscapeParams);
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

  return (
    <div className="position-relative">
      <div className="p-3 shadow-sm bg-white">
        <Card>
          <Header>
            <Nav variant="pills" defaultActiveKey="#mutationalProfiles">
              {[
                { title: 'Signature Exploring', id: 'signatureExploring' },
                { title: 'Exposure Exploring', id: 'exposureExploring' },
              ].map(({ title, id }) => {
                return (
                  <Item key={id}>
                    <Link
                      active={displayTab == id}
                      onClick={() => dispatchExploring({ displayTab: id })}
                    >
                      {title}
                    </Link>
                  </Item>
                );
              })}
            </Nav>
          </Header>
          <Body
            style={{
              display: displayTab == 'signatureExploring' ? 'block' : 'none',
            }}
          >
            <SignatureExploring />
          </Body>
          <Body
            style={{
              display: displayTab == 'exposureExploring' ? 'block' : 'none',
            }}
          >
            <ExposureExploring />
          </Body>
        </Card>
      </div>
    </div>
  );
}
