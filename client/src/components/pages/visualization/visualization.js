import React, { useEffect, useState } from 'react';
import { Form, Row, Col, Nav, Button } from 'react-bootstrap';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '../../controls/sidebar-container/sidebar-container';
import UserForm from './userForm';
import PublicForm from './publicForm';
import Instructions from '../visualization/instructions';
import ProfilerSummary from './profilerSummary';
import MutationalProfiles from './mutationalProfiles';
import TreeAndLeaf from './treeAndLeaf';
import CosineSimilarity from './cosineSimilarity';
import MutationalPattern from './mutationalPattern';
import ProfileComparison from './profileComparison';
import PCA from './pca';
import Kataegis from './kataegis';
import Download from './download';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../services/store/visualization';
import { actions as modalActions } from '../../../services/store/modal';
import {
  defaultProfile,
  defaultMatrix,
  defaultFilter,
} from '../../../services/utils';
import './visualization.scss';

const actions = { ...visualizationActions, ...modalActions };
const { Group, Label, Check } = Form;

export default function Visualization({ match }) {
  const dispatch = useDispatch();
  const visualization = useSelector((state) => state.visualization);
  const mergeState = (state) =>
    dispatch(actions.mergeVisualization({ state: state }));
  const mergeMutationalProfiles = (state) =>
    dispatch(actions.mergeVisualization({ mutationalProfiles: state }));
  const mergeKataegis = (state) =>
    dispatch(actions.mergeVisualization({ kataegis: state }));
  const mergeCosineSimilarity = (state) =>
    dispatch(actions.mergeVisualization({ cosineSimilarity: state }));
  const mergeProfileComparison = (state) =>
    dispatch(actions.mergeVisualization({ profileComparison: state }));
  const mergePCA = (state) =>
    dispatch(actions.mergeVisualization({ pca: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const {
    openSidebar,
    loading,
    source,
    submitted,
    queueExpired,
    error,
    projectID,
    displayTab,
    svgList,
    matrixList,
  } = visualization.state;
  const mutationalProfiles = visualization.mutationalProfiles;
  const { signatureSetOptions } = visualization.pca;

  const { type, id } = match.params;

  // when retrieving queued result, update id in store
  useEffect(() => {
    if (id && !loading.active && !submitted && !projectID) {
      if (type == 'queue') {
        loadQueueResult(id);
      } else if (type == 'example') {
        loadExample(id);
      }
    }
  }, [id, loading.active]);

  // get mapping of plots after retrieving projectID
  useEffect(() => {
    if (source == 'user') {
      if (projectID && !Object.keys(svgList).length) {
        getResults();
      } else if (Object.keys(svgList).length && !signatureSetOptions.length) {
        loadData();
      }
    } else {
      if (
        Object.keys(svgList).length > 0 &&
        !mutationalProfiles.filtered.length
      )
        mapPublicData();
    }
  }, [svgList, projectID]);

  // reload summary information
  async function getResults() {
    const response = await fetch(`api/getResults`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ projectID: projectID }),
    });

    if (response.ok) {
      const { svgList, statistics, matrixList, downloads } =
        await response.json();
      mergeState({
        svgList: svgList,
        statistics: statistics,
        matrixList: matrixList,
        downloads: downloads,
      });
    } else {
      mergeError(await response.json());
    }
  }

  // retrieve mapping of samples to plots from svgList file
  async function loadData() {
    mergeState({
      loading: {
        active: true,
        // content: 'Putting Data Into Session',
        // showIndicator: true,
      },
    });

    // Mutational Profiles
    const nameOptions = [...new Set(svgList.map((d) => d.Sample_Name))];
    const selectName = nameOptions[0];
    const filteredPlots = svgList.filter(
      (row) => row.Sample_Name == selectName
    );

    const filteredProfileOptions = [
      ...new Set(
        filteredPlots.map((row) => row.Profile_Type).sort((a, b) => a - b)
      ),
    ];

    const profile = defaultProfile(filteredProfileOptions);

    const filteredMatrixOptions = [
      ...new Set(
        filteredPlots
          .filter((row) => row.Profile_Type == profile)
          .map((row) => row.Matrix_Size)
      ),
    ].sort((a, b) => a - b);

    const matrix = defaultMatrix(profile, filteredMatrixOptions);

    const filteredFilterOptions = [
      ...new Set(
        filteredPlots
          .filter(
            (row) => row.Profile_Type == profile && row.Matrix_Size == matrix
          )
          .map((row) => row.Filter)
          .sort((a, b) => a - b)
      ),
    ];

    const filter = defaultFilter(filteredFilterOptions);

    const filteredMatrixList = [
      ...new Set(
        matrixList
          .filter((row) => row.Profile_Type == profile)
          .map((row) => row.Matrix_Size)
          .sort((a, b) => a - b)
      ),
    ];

    mergeMutationalProfiles({
      filtered: filteredPlots,
      nameOptions: nameOptions,
      profileOptions: filteredProfileOptions,
      matrixOptions: filteredMatrixOptions,
      filterOptions: filteredFilterOptions,
      selectName: selectName,
      selectProfile: profile,
      selectMatrix: matrix,
      selectFilter: filter,
    });

    // Cosine Similarity - Profile Comparison - PCA - Kataegis
    const sampleNameOptions = [
      ...new Set(
        svgList.map((row) => {
          if (row.Filter != 'NA') return `${row.Sample_Name}@${row.Filter}`;
          else return row.Sample_Name;
        })
      ),
    ];
    const profileOptions = [...new Set(svgList.map((row) => row.Profile_Type))];

    const selectProfile = defaultProfile(profileOptions);
    const selectMatrix = defaultMatrix(selectProfile, filteredMatrixOptions);

    const { output: refSignatureSetOptions } = await (
      await getRefSigOptions(selectProfile)
    ).json();

    mergeCosineSimilarity({
      withinProfileType: selectProfile,
      refProfileType: selectProfile,
      refSignatureSet: refSignatureSetOptions[0],
      refSignatureSetOptions: refSignatureSetOptions,
      withinMatrixSize: selectMatrix,
      withinMatrixOptions: filteredMatrixList,
      userProfileType: selectProfile,
      userMatrixSize: selectMatrix,
      userMatrixOptions: filteredMatrixOptions,
    });

    mergeProfileComparison({
      withinProfileType: selectProfile,
      withinSampleName1: sampleNameOptions[0],
      withinSampleName2: sampleNameOptions[1],
      sampleOptions: sampleNameOptions,
      refProfileType: selectProfile,
      refSampleName: sampleNameOptions[0],
      refSignatureSet: refSignatureSetOptions[0],
      refSignatureSetOptions: refSignatureSetOptions,
      userProfileType: selectProfile,
      userMatrixSize: selectMatrix,
      userMatrixOptions: filteredMatrixOptions,
      userSampleName: sampleNameOptions[0],
    });

    mergePCA({
      profileType: selectProfile,
      signatureSet: refSignatureSetOptions[0],
      signatureSetOptions: refSignatureSetOptions,
      userProfileType: selectProfile,
      userMatrixSize: selectMatrix,
      userMatrixOptions: filteredMatrixOptions,
    });

    mergeKataegis({
      sample: sampleNameOptions[0],
      sampleOptions: sampleNameOptions,
    });

    mergeState({
      profileOptions: profileOptions,
      loading: {
        active: false,
      },
      openSidebar: false,
    });
  }

  // retrieve mapping of samples to plots from svgList file
  async function mapPublicData() {
    mergeState({
      loading: {
        active: true,
        // content: 'Putting Public Data Into Session',
        // showIndicator: true,
      },
    });

    // Mutational Profiles
    const nameOptions = [...new Set(svgList.map((row) => row.Sample))];
    const selectName = mutationalProfiles.selectName || nameOptions[0];
    const profileOptions = [
      ...new Set(svgList.map((row) => row.Profile.match(/[a-z]+/gi)[0])),
    ];
    const profile = defaultProfile(profileOptions);

    const filteredMatrixOptions = [
      ...new Set(
        svgList
          .filter(
            (row) =>
              row.Sample == selectName && row.Profile.indexOf(profile) > -1
          )
          .map(({ Profile }) => Profile.match(/\d+/gi)[0])
      ),
    ].sort((a, b) => a - b);

    mergeState({ profileOptions: profileOptions });
    // Cosine Similarity - Profile Comparison - PCA
    const selectProfile = defaultProfile(profileOptions);
    const selectMatrix = defaultMatrix(selectProfile, filteredMatrixOptions);

    const { output: refSignatureSetOptions } = await (
      await getRefSigOptions(selectProfile)
    ).json();

    mergeCosineSimilarity({
      withinProfileType: selectProfile,
      refProfileType: selectProfile,
      refSignatureSet: refSignatureSetOptions[0],
      refSignatureSetOptions: refSignatureSetOptions,
      withinMatrixSize: selectMatrix,
      withinMatrixOptions: filteredMatrixOptions,
    });

    mergeProfileComparison({
      withinProfileType: selectProfile,
      withinSampleName1: nameOptions[0],
      withinSampleName2: nameOptions[1],
      sampleOptions: nameOptions,
      refProfileType: selectProfile,
      refSampleName: nameOptions[0],
      refSignatureSet: refSignatureSetOptions[0],
      refSignatureSetOptions: refSignatureSetOptions,
    });

    mergePCA({
      profileType: selectProfile,
      signatureSet: refSignatureSetOptions[0],
      signatureSetOptions: refSignatureSetOptions,
    });

    mergeState({
      loading: {
        active: false,
      },
      submitted: true,
      displayTab: 'profilerSummary',
      openSidebar: false,
    });
  }

  function submitR(fn, args) {
    return fetch(`api/visualizationWrapper`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ fn: fn, args: args, projectID: projectID }),
    });
  }

  function getRefSigOptions(profileType) {
    return fetch(`api/visualizationWrapper`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        fn: 'getReferenceSignatureSets',
        args: { profileType: profileType },
      }),
    });
  }

  async function loadQueueResult(id) {
    mergeState({
      loading: {
        active: true,
        // content: 'Loading Queued Result',
        // showIndicator: true,
      },
    });
    try {
      const { args, visualization, timestamp } = await (
        await fetch(`api/getQueueResults/${id}`)
      ).json();
      dispatch(actions.mergeVisualization(visualization));
    } catch (error) {
      mergeState({
        queueExpired: true,
      });
    }
    mergeState({
      loading: { active: false },
      submitted: true,
      displayTab: 'profilerSummary',
      openSidebar: false,
    });
  }

  async function loadExample(id) {
    mergeState({
      loading: {
        active: true,
        // content: 'Loading Example',
        // showIndicator: true,
      },
    });
    try {
      const { projectID, state: visualizationStore } = await (
        await fetch(`api/getVisExample/${id}`)
      ).json();
      dispatch(
        actions.mergeVisualization({
          ...visualizationStore,
          state: { ...visualizationStore.state, projectID: projectID },
        })
      );
    } catch (error) {
      mergeError(`Example: ${id} does not exist`);
    }
    mergeState({
      loading: { active: false },
      submitted: true,
      displayTab: 'profilerSummary',
      openSidebar: false,
    });
  }

  const tabs = [
    {
      name: 'Instructions',
      id: 'instructions',
      component: <Instructions />,
    },
    {
      name: 'Profiler Summary',
      id: 'profilerSummary',
      component: <ProfilerSummary submitR={(fn, args) => submitR(fn, args)} />,
    },
    {
      name: 'Mutational Profiles',
      id: 'mutationalProfiles',
      component: <MutationalProfiles />,
    },
    {
      name: 'Tree And Leaf',
      id: 'treeAndLeaf',
      component: <TreeAndLeaf />,
    },
    {
      name: 'Cosine Similarity',
      id: 'cosineSimilarity',
      component: (
        <CosineSimilarity
          getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
          submitR={(fn, args) => submitR(fn, args)}
        />
      ),
    },
    {
      name: 'Mutational Pattern Enrichment Analysis',
      id: 'mutationalPattern',
      component: (
        <MutationalPattern submitR={(fn, args) => submitR(fn, args)} />
      ),
    },
    {
      name: 'Profile Comparison',
      id: 'profileComparison',
      component: (
        <ProfileComparison
          getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
          submitR={(fn, args) => submitR(fn, args)}
        />
      ),
    },
    {
      name: 'PCA',
      id: 'pca',
      component: (
        <PCA
          getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
          submitR={(fn, args) => submitR(fn, args)}
        />
      ),
    },
    {
      name: 'Kataegis Identification',
      id: 'kataegisIdentification',
      component: <Kataegis submitR={(fn, args) => submitR(fn, args)} />,
    },
    {
      name: 'Download',
      id: 'download',
      component: <Download />,
    },
  ];

  return (
    <div className="position-relative">
      <div className="mx-3">
        <div className="mx-3 bg-white border border-top-0">
          {/* for desktops and tablets */}
          <div className="d-none d-md-block">
            <Nav defaultActiveKey="profilerSummary">
              {tabs.map(({ name, id }) => (
                <div key={id} className="d-inline-block">
                  <Button
                    variant="link"
                    className={`secondary-navlinks px-3 py-1 d-inline-block border-0 ${
                      id == displayTab ? 'active-secondary-navlinks' : ''
                    }`}
                    active={id == displayTab && submitted}
                    disabled={id != 'instructions' && !submitted}
                    style={{
                      textDecoration: 'none',
                      fontSize: '12pt',
                      color: 'black',
                      fontWeight: '500',
                    }}
                    onClick={() => mergeState({ displayTab: id })}
                  >
                    {name}
                  </Button>
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
                  <Button
                    variant="link"
                    className={
                      id == displayTab && Object.keys(svgList).length
                        ? 'secondary-navlinks px-3 py-1 d-inline-block border-0 active-secondary-navlinks'
                        : 'secondary-navlinks px-3 py-1 d-inline-block border-0'
                    }
                    style={{
                      textDecoration: 'none',
                      fontSize: '12pt',
                      color: 'black',
                      fontWeight: '500',
                    }}
                    onClick={() => mergeState({ displayTab: id })}
                  >
                    {name}
                  </Button>
                  <div className="d-md-none w-100"></div>
                </div>
              ))}
            </Nav>
          </div>
        </div>
      </div>
      <SidebarContainer
        className="m-3"
        collapsed={!openSidebar}
        onCollapsed={() => mergeState({ openSidebar: !openSidebar })}
      >
        <SidebarPanel className="col-lg-3 col-md-5">
          <div className="p-3 bg-white border rounded">
            <Row>
              <Col lg="auto">
                <Group>
                  <Label className="mr-4">Data Source</Label>
                  <Check inline id="radioPublic">
                    <Check.Input
                      disabled={submitted || loading.active}
                      type="radio"
                      value="public"
                      checked={source == 'public'}
                      onChange={(e) => mergeState({ source: 'public' })}
                    />
                    <Check.Label className="font-weight-normal">
                      Public
                    </Check.Label>
                  </Check>
                  <Check inline id="radioUser">
                    <Check.Input
                      disabled={submitted || loading.active}
                      type="radio"
                      value="user"
                      checked={source == 'user'}
                      onChange={(e) => mergeState({ source: 'user' })}
                    />
                    <Check.Label className="font-weight-normal">
                      User
                    </Check.Label>
                  </Check>
                </Group>
              </Col>
            </Row>
            <Row
              style={{
                display: source == 'user' ? 'block' : 'none',
              }}
            >
              <Col lg="auto" className="w-100">
                <UserForm />
              </Col>
            </Row>
            <Row style={{ display: source == 'public' ? 'block' : 'none' }}>
              <Col lg="auto" className="w-100">
                <PublicForm />
              </Col>
            </Row>
          </div>
          <hr className="d-lg-none" style={{ opacity: 0 }}></hr>
        </SidebarPanel>
        <MainPanel className="col-lg-9 col-md-7">
          {queueExpired && (
            <div className="bg-warning mb-3 p-3 border rounded">
              <span>Queue results have expired</span>
            </div>
          )}
          <div>
            <LoadingOverlay
              active={loading.active}
              content={loading.content}
              showIndicator={loading.showIndicator}
            />
            <div style={{ minHeight: '500px' }}>
              {tabs.filter((tab) => tab.id == displayTab)[0].component}
            </div>
          </div>
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
