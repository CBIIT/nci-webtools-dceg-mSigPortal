import React, { useEffect } from 'react';
import axios from 'axios';
import { Form, Row, Col, Nav, Button } from 'react-bootstrap';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '../../controls/sidebar-container/sidebar-container';
import UserForm from './userForm/userForm';
import PublicForm from './publicForm/publicForm';
import Instructions from '../visualization/instructions';
import ProfilerSummary from './profilerSummary/profilerSummary';
import MutationalProfiles from './mutationalProfiles';
import TreeAndLeaf from './treeLeaf/treeLeaf';
import MutationalProfiles2 from './mutationalProfiles/mutProfiles';
import CosineSimilarity from './cosineSimilarity';
import CosineSimilarity2 from './cosineSimilarity/cosineSimilarity';
import MutationalPattern from './mutationalPattern/mutationalPattern';
import ProfileComparison from './profileComparison';
import ProfileComparison2 from './profileComparison/profileComparison';
import PCA from './pca';
import PCA2 from './pca/pca';
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

import MultationalProfilesTest from './test/multationalProfilesTest';

const actions = { ...visualizationActions, ...modalActions };
const { Group, Label, Check } = Form;

export default function Visualization({ match }) {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.visualization);
  const mergeState = (state) =>
    dispatch(actions.mergeVisualization({ main: state }));
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
    samples,
  } = store.main;
  const mutationalProfiles = store.mutationalProfiles;
  const { signatureSetOptions } = store.pca;

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
    }
  }, [svgList, projectID]);

  // switch to first tab after fetching samples
  useEffect(() => {
    if (samples.length && projectID) {
      mergeState({ displayTab: 'profilerSummary' });
    }
  }, [samples, projectID]);

  // reload summary information
  async function getResults() {
    const response = await fetch(`web/getResults`, {
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
    const nameOptions = [...new Set(svgList.map((d) => d.sample))];
    const selectName = nameOptions[0];
    const filteredPlots = svgList.filter((e) => e.sample == selectName);

    const filteredProfileOptions = [
      ...new Set(filteredPlots.map((e) => e.profileType).sort((a, b) => a - b)),
    ];

    const profile = defaultProfile(filteredProfileOptions);

    const filteredMatrixOptions = [
      ...new Set(
        filteredPlots
          .filter((e) => e.profileType == profile)
          .map((e) => e.matrixSize)
      ),
    ].sort((a, b) => a - b);

    const matrix = defaultMatrix(profile, filteredMatrixOptions);

    const filteredFilterOptions = [
      ...new Set(
        filteredPlots
          .filter((e) => e.profileType == profile && e.matrixSize == matrix)
          .map((e) => e.Filter)
          .sort((a, b) => a - b)
      ),
    ];

    const filter = defaultFilter(filteredFilterOptions);

    const filteredMatrixList = [
      ...new Set(
        matrixList
          .filter((e) => e.profileType == profile)
          .map((e) => e.matrixSize)
          .sort((a, b) => a - b)
      ),
    ];

    // mergeMutationalProfiles({
    //   filtered: filteredPlots,
    //   nameOptions: nameOptions,
    //   profileOptions: filteredProfileOptions,
    //   matrixOptions: filteredMatrixOptions,
    //   filterOptions: filteredFilterOptions,
    //   selectName: selectName,
    //   selectProfile: profile,
    //   selectMatrix: matrix,
    //   selectFilter: filter,
    // });

    // Cosine Similarity - Profile Comparison - PCA - Kataegis
    const sampleNameOptions = [
      ...new Set(
        svgList.map((e) => {
          if (e.Filter != 'NA') return `${e.sample}@${e.Filter}`;
          else return e.samples;
        })
      ),
    ];
    const profileOptions = [...new Set(svgList.map((e) => e.profileType))];

    const selectProfile = defaultProfile(profileOptions);
    const selectMatrix = defaultMatrix(selectProfile, filteredMatrixOptions);

    // const { output: refSignatureSetOptions } = await (
    //   await getRefSigOptions(selectProfile)
    // ).json();

    mergeCosineSimilarity({
      withinProfileType: selectProfile,
      refProfileType: selectProfile,
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
      // refSignatureSet: refSignatureSetOptions[0],
      // refSignatureSetOptions: refSignatureSetOptions,
      userProfileType: selectProfile,
      userMatrixSize: selectMatrix,
      userMatrixOptions: filteredMatrixOptions,
      userSampleName: sampleNameOptions[0],
    });

    mergePCA({
      profileType: selectProfile,
      // signatureSet: refSignatureSetOptions[0],
      // signatureSetOptions: refSignatureSetOptions,
      userProfileType: selectProfile,
      userMatrixSize: selectMatrix,
      userMatrixOptions: filteredMatrixOptions,
    });

    mergeKataegis({
      sample: sampleNameOptions[0],
      sampleOptions: sampleNameOptions,
    });

    mergeState({
      displayTab: 'profilerSummary',
      profileOptions: profileOptions,
      loading: {
        active: false,
      },
      openSidebar: false,
    });
  }

  function submitR(fn, args) {
    return fetch(`web/visualizationWrapper`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ fn: fn, args: args, projectID: projectID }),
    });
  }

  // function getRefSigOptions(profileType) {
  //   return fetch(`web/visualizationWrapper`, {
  //     method: 'POST',
  //     headers: {
  //       Accept: 'application/json',
  //       'Content-Type': 'application/json',
  //     },
  //     body: JSON.stringify({
  //       fn: 'getReferenceSignatureSets',
  //       args: { profileType: profileType },
  //     }),
  //   });
  // }

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
        await fetch(`web/getQueueResults/${id}`)
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
        await fetch(`web/getVisExample/${id}`)
      ).json();
      dispatch(
        actions.mergeVisualization({
          ...visualizationStore,
          state: { ...visualizationStore.main, projectID: projectID },
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
      component: <ProfilerSummary />,
    },
    {
      name: 'Mutational Profiles',
      id: 'mutationalProfiles',
      component:
        source == 'user' ? <MutationalProfiles /> : <MutationalProfiles2 />,
    },
    {
      name: 'Tree and Leaf',
      id: 'treeAndLeaf',
      component: <TreeAndLeaf />,
    },
    {
      name: 'Cosine Similarity',
      id: 'cosineSimilarity',
      component:
        source == 'user' ? (
          <CosineSimilarity
            // getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
            submitR={(fn, args) => submitR(fn, args)}
          />
        ) : (
          <CosineSimilarity2 />
        ),
    },
    {
      name: 'Mutational Pattern Enrichment Analysis',
      id: 'mutationalPattern',
      component: <MutationalPattern />,
    },
    {
      name: 'Profile Comparison',
      id: 'profileComparison',
      component:
        source == 'user' ? (
          <ProfileComparison
            // getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
            submitR={(fn, args) => submitR(fn, args)}
          />
        ) : (
          <ProfileComparison2 />
        ),
    },

    {
      name: 'PCA',
      id: 'pca',
      component:
        source == 'user' ? (
          <PCA
            // getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
            submitR={(fn, args) => submitR(fn, args)}
          />
        ) : (
          <PCA2 />
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
    {
      name: 'Multational Profiles Test',
      id: 'multationalProfileTest',
      component: <MultationalProfilesTest />,
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
          <div className="e d-md-none">
            <Nav defaultActiveKey="summary">
              {tabs.map(({ name, id }) => (
                <div key={id} className="col-12 text-center">
                  <Button
                    variant="link"
                    className={
                      id == displayTab &&
                      (samples.length || Object.keys(svgList).length)
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
          {error && (
            <div className="bg-danger text-white mb-3 p-3 border rounded">
              {error}
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
              {/* {tabs.map((tab) => (
                <div className={tab.id == displayTab ? 'd-block' : 'd-none'}>
                  {tab.component}
                </div>
              ))} */}
            </div>
          </div>
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
