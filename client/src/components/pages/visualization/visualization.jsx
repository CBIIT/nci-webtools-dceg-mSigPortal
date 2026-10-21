import { useCallback, useEffect, useMemo } from 'react';
import { Form, Row, Col, Nav, Button, Alert } from 'react-bootstrap';
import { useParams } from 'react-router-dom';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '@/components/controls/sidebar-container/sidebar-container';
import UserForm from './userForm/userForm';
import PublicForm from './publicForm/publicForm';
import Instructions from '@/components/pages/visualization/instructions';
import ProfilerSummary from './profilerSummary/profilerSummary';
import MutationalProfiles from './mutationalProfiles/mutProfiles';
import TreeAndLeaf from './treeLeaf/treeLeaf';
import CosineSimilarity from './cosineSimilarity/cosineSimilarity';
import MutationalPattern from './mutationalPattern/mutationalPattern';
import ProfileComparison from './profileComparison/profileComparison';
import PCA from './pca/pca';
import ClusteredIdentification from './clustered/clustered';
import Download from './download';
import { LoadingOverlay } from '@/components/controls/loading-overlay/loading-overlay';
import TabPanel from '@/components/controls/tab/tab';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '@/services/store/visualization';
import { actions as modalActions } from '@/services/store/modal';
import './visualization.scss';
import { useRefreshQuery } from './userForm/apiSlice';

const actions = { ...visualizationActions, ...modalActions };
const { Group, Label, Check } = Form;

export default function Visualization() {
  const dispatch = useDispatch();
  const mergeState = (state) =>
    dispatch(actions.mergeVisualization({ main: state }));
  const { main, publicForm, mutationalProfiles } = useSelector(
    (state) => state.visualization
  );
  const {
    openSidebar,
    loading,
    source,
    submitted,
    error,
    displayTab,
    ...state
  } = main;

  const id = useParams().id || state.id;
  const {
    data: jobInfo,
    error: refreshError,
    refetch: refreshJob,
  } = useRefreshQuery(id, {
    skip: !id,
  });
  const status = jobInfo?.status;
  const params = jobInfo?.params;
  const isDone = ['COMPLETED', 'FAILED'].includes(status?.status);
  const completed =
    (source === 'user' && status?.status === 'COMPLETED') ||
    (source === 'public' && submitted);

  const tabState = useMemo(
    () => ({ ...publicForm, ...main, id, source }),
    [publicForm, main, id, source]
  );
  const mutationalProfilesState = useMemo(
    () => ({ ...tabState, mutationalProfiles }),
    [tabState, mutationalProfiles]
  );
  const clusterState = useMemo(
    () => ({ ...publicForm, ...main, params, id, source }),
    [publicForm, main, params, id, source]
  );

  const refreshState = useCallback(() => {
    refreshJob();
  }, [refreshJob]);

  // parse py script error message and return a user friendly message
  function errorMessage() {
    const error = status?.error;
    if (error) {
      if (error.includes('Error 2727')) {
        return error;
      } else {
        return 'An error has occurred. Please review your inputs and try again.';
      }
    } else {
      return 'INTERNAL ERROR';
    }
  }

  // refresh job status every minute
  useEffect(() => {
    const interval = setInterval(refreshState, 1000 * 60);
    if (isDone || refreshError) clearInterval(interval);
    return () => clearInterval(interval);
  }, [isDone, refreshError, refreshState]);
  // switch to first tab when job is complete
  useEffect(() => {
    if (status && status.status === 'COMPLETED' && displayTab == 'instructions')
      mergeState({
        displayTab: 'profilerSummary',
        openSidebar: false,
        source: 'user',
      });
  }, [status]);

  const tabs = [
    {
      name: 'Instructions',
      id: 'instructions',
      disabled: false,
    },
    {
      name: 'Profiler Summary',
      id: 'profilerSummary',
      disabled: !completed,
    },
    {
      name: 'Mutational Profiles',
      id: 'mutationalProfiles',
      disabled: !completed,
    },
    {
      name: 'Tree and Leaf',
      id: 'treeAndLeaf',
      disabled: !completed,
    },
    {
      name: 'Cosine Similarity',
      id: 'cosineSimilarity',
      disabled: !completed,
    },
    {
      name: 'Profile Comparison',
      id: 'profileComparison',
      disabled: !completed,
    },
    {
      name: 'PCA',
      id: 'pca',
      disabled: !completed,
    },
    {
      name: 'Mutational Pattern Enrichment Analysis',
      id: 'mutationalPattern',
      disabled: !completed,
    },

    source == 'user' && {
      name: 'Clustered Mutations Identification',
      id: 'cluster',
      disabled: !completed,
    },

    // {
    //   name: 'Download',
    //   id: 'download',
    //   component: <Download />,
    // },
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
                .map(({ name, id, disabled }, i) => (
                  <div key={name} className="d-inline-block">
                    <Button
                      variant="link"
                      className={`secondary-navlinks px-3 py-1 d-inline-block border-0 text-exploration rounded-0 ${
                        id === displayTab
                          ? 'bg-visualization text-white'
                          : 'text-visualization'
                      }`}
                      disabled={disabled}
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
                  </div>
                ))}
            </Nav>
          </div>

          {/* for mobile devices */}
          <div className="e d-md-none">
            <Nav defaultActiveKey="summary">
              {tabs
                .filter((e) => e)
                .map(({ name, id, disabled }, i) => {
                  if (name)
                    return (
                      <div key={name} className="col-12 text-center">
                        <Button
                          variant="link"
                          className={
                            id === displayTab
                              ? 'secondary-navlinks px-3 py-1 d-inline-block border-0 bg-visualization text-white rounded-0'
                              : 'secondary-navlinks px-3 py-1 d-inline-block border-0 rounded-0'
                          }
                          disabled={disabled}
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
                    );
                })}
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
          {refreshError && <Alert variant="danger">Results expired</Alert>}
          {status && status.status === 'SUBMITTED' && (
            <div className="border rounded bg-white mb-3 p-3">
              <p>
                Your job has been submitted. You will receive an email once it
                is complete.
              </p>
            </div>
          )}
          {status && status.status === 'IN_PROGRESS' && (
            <div className="border rounded bg-white mb-3 p-3">
              <p>Your analysis is currently in progress.</p>
              <LoadingOverlay active={true} />
            </div>
          )}
          {status && status.status === 'FAILED' && (
            <div className="border rounded bg-white mb-3 p-3">
              <p>
                An error occurred during calculation:
                <br />
                {errorMessage()}
                {source == 'user' && (
                  <div>
                    <Button
                      download
                      variant="link"
                      className="p-0"
                      href={`api/data/output/${id}/profiler_extraction_log.txt`}
                    >
                      Download Job Log
                    </Button>
                  </div>
                )}
                <br />
                Please contact the site administrator for assistance if this
                issue persists.
              </p>
            </div>
          )}

          <TabPanel active={displayTab === 'instructions'}>
            <Instructions />
          </TabPanel>
          {completed && (
            <>
              <TabPanel active={displayTab === 'profilerSummary'}>
                <ProfilerSummary state={tabState} />
              </TabPanel>
              <TabPanel active={displayTab === 'mutationalProfiles'}>
                <MutationalProfiles state={mutationalProfilesState} />
              </TabPanel>
              <TabPanel active={displayTab === 'treeAndLeaf'}>
                <TreeAndLeaf state={tabState} />
              </TabPanel>
              <TabPanel active={displayTab === 'cosineSimilarity'}>
                <CosineSimilarity state={tabState} />
              </TabPanel>
              <TabPanel active={displayTab === 'profileComparison'}>
                <ProfileComparison state={tabState} />
              </TabPanel>
              <TabPanel active={displayTab === 'pca'}>
                <PCA state={tabState} />
              </TabPanel>
              <TabPanel active={displayTab === 'mutationalPattern'}>
                <MutationalPattern state={tabState} />
              </TabPanel>
              {source == 'user' && (
                <TabPanel active={displayTab === 'cluster'}>
                  <ClusteredIdentification state={clusterState} />
                </TabPanel>
              )}
              <TabPanel active={displayTab === 'download'}>
                <Download state={tabState} />
              </TabPanel>
            </>
          )}
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
