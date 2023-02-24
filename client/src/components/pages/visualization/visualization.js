import React, { useEffect } from 'react';
import { Form, Row, Col, Nav, Button } from 'react-bootstrap';
import { useParams } from 'react-router-dom';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '../../controls/sidebar-container/sidebar-container';
import UserForm from './userForm/userForm';
import PublicForm from './publicForm/publicForm';
import Instructions from '../visualization/instructions';
import ProfilerSummary from './profilerSummary/profilerSummary';
import MutationalProfiles from './mutationalProfiles/mutProfiles';
import TreeAndLeaf from './treeLeaf/treeLeaf';
import CosineSimilarity from './cosineSimilarity/cosineSimilarity';
import MutationalPattern from './mutationalPattern/mutationalPattern';
import ProfileComparison from './profileComparison/profileComparison';
import PCA from './pca/pca';
import ClusteredIdentification from './clustered/clustered';
import Download from './download';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../services/store/visualization';
import { actions as modalActions } from '../../../services/store/modal';
import './visualization.scss';

const actions = { ...visualizationActions, ...modalActions };
const { Group, Label, Check } = Form;

export default function Visualization() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.visualization);
  const mergeState = (state) =>
    dispatch(actions.mergeVisualization({ main: state }));

  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const {
    openSidebar,
    loading,
    source,
    submitted,
    queueExpired,
    error,
    displayTab,
    svgList,
    matrixList,
    matrixData,
  } = store.main;

  const urlParams = useParams();
  const { type } = urlParams;
  const id = store.main.id || urlParams.id;

  // when retrieving queued result, update id in store
  useEffect(() => {
    if (id && !loading.active && !submitted && !id) {
      if (type == 'queue') {
        loadQueueResult(id);
      } else if (type == 'example') {
        loadExample(id);
      }
    }
  }, [id, loading.active]);

  // switch to first tab after fetching samples
  useEffect(() => {
    if (matrixData.length) {
      mergeState({ displayTab: 'profilerSummary', openSidebar: false });
    }
  }, [matrixData]);

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
      // displayTab: 'profilerSummary',
      // openSidebar: false,
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
      const { id, state: visualizationStore } = await (
        await fetch(`web/getVisExample/${id}`)
      ).json();
      dispatch(
        actions.mergeVisualization({
          ...visualizationStore,
          state: { ...visualizationStore.main, id },
        })
      );
    } catch (error) {
      mergeError(`Example: ${id} does not exist`);
    }
    mergeState({
      loading: { active: false },
      submitted: true,
      // displayTab: 'profilerSummary',
      // openSidebar: false,
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
      component: <MutationalProfiles />,
    },
    {
      name: 'Tree and Leaf',
      id: 'treeAndLeaf',
      component: <TreeAndLeaf />,
    },
    {
      name: 'Cosine Similarity',
      id: 'cosineSimilarity',
      component: <CosineSimilarity />,
    },
    {
      name: 'Profile Comparison',
      id: 'profileComparison',
      component: <ProfileComparison />,
    },
    {
      name: 'PCA',
      id: 'pca',
      component: <PCA />,
    },
    {
      name: 'Mutational Pattern Enrichment Analysis',
      id: 'mutationalPattern',
      component: <MutationalPattern />,
    },

    source == 'user' && {
      name: 'Clustered Mutations Identification',
      id: 'cluster',
      component: <ClusteredIdentification />,
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
              {tabs
                .filter((e) => e)
                .map(({ name, id }) => (
                  <div key={id} className="d-inline-block">
                    <Button
                      variant="link"
                      className={`secondary-navlinks px-3 py-1 d-inline-block border-0 rounded-0 ${
                        id == displayTab
                          ? 'bg-visualization text-white'
                          : 'text-visualization'
                      }`}
                      active={id == displayTab}
                      disabled={
                        id != 'instructions' &&
                        !(source == 'public'
                          ? matrixData.length
                          : matrixList.length)
                      }
                      style={{
                        textDecoration: 'none',
                        fontSize: '12pt',
                        color: '#3a7867',
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
                <div key={id + 'm'} className="col-12 text-center">
                  <Button
                    variant="link"
                    className={
                      id == displayTab &&
                      (matrixData.length || Object.keys(svgList).length)
                        ? 'secondary-navlinks px-3 py-1 d-inline-block border-0 bg-visualization text-white rounded-0'
                        : 'secondary-navlinks px-3 py-1 d-inline-block border-0 rounded-0'
                    }
                    style={{
                      textDecoration: 'none',
                      fontSize: '12pt',
                      color: '#3a7867',
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
