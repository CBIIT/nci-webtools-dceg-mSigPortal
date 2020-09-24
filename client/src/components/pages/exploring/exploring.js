import React, { useEffect } from 'react';
import { Card, Nav, Row, Col, Form } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '../../controls/sidebar-container/sidebar-container';
import SignatureExploring from './signature/signatureExploring';
import ExposureExploring from './exposure/exposureExploring';
import {
  dispatchError,
  dispatchExploring,
  dispatchExpMutationalProfiles,
  dispatchExpCosineSimilarity,
  dispatchExpMutationalSigComparison,
  dispatchExpTumor,
  dispatchExpActivity,
  dispatchExpDecomposition,
  dispatchExpLandscape,
  dispatchExpPrevalence,
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
    source,
    openSidebar,
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
    dispatchExpTumor({ loading: true });
    dispatchExpActivity({ loading: true });
    dispatchExpDecomposition({ loading: true });
    dispatchExpLandscape({ loading: true });
    dispatchExpPrevalence({ loading: true });

    try {
      let [publicData, referenceSignatures] = await Promise.all([
        fetch(`${rootURL}getPublicDataOptions`, {
          method: 'POST',
          headers: {
            Accept: 'application/json',
            'Content-Type': 'application/json',
          },
        }),
        getReferenceSignatureData(['Signature_set_name']),
      ]);

      if (!publicData.ok) {
        // console.log(await response.json());
      } else {
        const data = await publicData.json();
        const studyOptions = [...new Set(data.map((data) => data.Study))];
        const strategyOptions = [
          ...new Set(
            data
              .filter((data) => data.Study == studyOptions[0])
              .map((data) => data.Dataset)
          ),
        ];
        const refSignatureSetOptions = [
          ...new Set(
            referenceSignatures.output.data.map((row) => row.Signature_set_name)
          ),
        ];

        const params = {
          study: studyOptions[0],
          studyOptions: studyOptions,
          strategy: strategyOptions[0],
          strategyOptions: strategyOptions,
          refSignatureSet: refSignatureSetOptions[0],
          refSignatureSetOptions: refSignatureSetOptions,
        };
        dispatchExploring({ publicDataOptions: data });
        dispatchExpTumor({ ...params });
        dispatchExpActivity({ ...params });
        dispatchExpDecomposition({ ...params });
        dispatchExpLandscape({ ...params });
        dispatchExpPrevalence({ ...params });
      }
    } catch (err) {
      dispatchError(err);
    }
    dispatchExpTumor({ loading: false });
    dispatchExpActivity({ loading: false });
    dispatchExpDecomposition({ loading: false });
    dispatchExpLandscape({ loading: false });
    dispatchExpPrevalence({ loading: false });
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
      <SidebarContainer
        className="m-3"
        collapsed={!openSidebar}
        onCollapsed={(e) => dispatchExploring({ openSidebar: !e })}
      >
        <SidebarPanel>
          <div className="p-3 shadow-sm bg-white">
            <Row>
              <Col sm="auto">
                <Group className="d-flex">
                  <Label className="mr-auto">
                    <h3 className="mb-2">Data Source</h3>
                  </Label>
                  <Check inline id="radioUser" className="ml-4">
                    <Check.Input
                      disabled={submitted}
                      type="radio"
                      value="user"
                      checked={source == 'user'}
                      onChange={(e) => dispatchExploring({ source: 'user' })}
                    />
                    <Check.Label className="font-weight-normal">
                      User
                    </Check.Label>
                  </Check>
                  <Check inline id="radioPublic">
                    <Check.Input
                      disabled={submitted}
                      type="radio"
                      value="public"
                      checked={source == 'public'}
                      onChange={(e) => dispatchExploring({ source: 'public' })}
                    />
                    <Check.Label className="font-weight-normal">
                      Public
                    </Check.Label>
                  </Check>
                </Group>
              </Col>
            </Row>
            <Row style={{ display: source == 'user' ? 'block' : 'none' }}>
              <Col sm="auto" className="w-100">
                {/* <UploadForm /> */}
              </Col>
            </Row>
            <Row style={{ display: source == 'public' ? 'block' : 'none' }}>
              <Col sm="auto" className="w-100">
                {/* <PublicForm /> */}
              </Col>
            </Row>
          </div>
          <hr className="d-lg-none" style={{ opacity: 0 }}></hr>
        </SidebarPanel>
        <MainPanel>
          <div
            className="p-3 shadow-sm bg-white"
            style={{ minHeight: '420px' }}
          >
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
                  display:
                    displayTab == 'signatureExploring' ? 'block' : 'none',
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
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
