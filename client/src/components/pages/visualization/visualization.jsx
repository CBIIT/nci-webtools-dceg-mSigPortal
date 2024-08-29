import { useCallback, useEffect } from 'react';
import { Form, Row, Col, Nav, Button, Alert } from 'react-bootstrap';
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

          <div className={displayTab === 'instructions' ? 'd-block' : 'd-none'}>
            <Instructions />
          </div>
          {completed && (
            <>
              <div
                className={
                  displayTab === 'profilerSummary' ? 'd-block' : 'd-none'
                }
              >
                <ProfilerSummary
                  state={{ ...publicForm, ...main, id, source }}
                />
              </div>
              <div
                className={
                  displayTab === 'mutationalProfiles' ? 'd-block' : 'd-none'
                }
              >
                <MutationalProfiles
                  state={{
                    ...publicForm,
                    ...main,
                    mutationalProfiles,
                    id,
                    source,
                  }}
                />
              </div>
              {displayTab === 'treeAndLeaf' && (
                source === 'public' 
                  ? <TreeAndLeaf state={{ ...publicForm, ...main, id, source }} /> 
                  : <div className="container-fluid bg-white border rounded p-3 text-center">
                      <Alert variant="warning">
                        The Tree and Leaf plot does not currently support user-provided data.
                      </Alert>
                    </div>
              )}
              <div
                className={
                  displayTab === 'cosineSimilarity' ? 'd-block' : 'd-none'
                }
              >
                <CosineSimilarity
                  state={{ ...publicForm, ...main, id, source }}
                />
              </div>
              <div
                className={
                  displayTab === 'profileComparison' ? 'd-block' : 'd-none'
                }
              >
                <ProfileComparison
                  state={{ ...publicForm, ...main, id, source }}
                />
              </div>
              <div className={displayTab === 'pca' ? 'd-block' : 'd-none'}>
                <PCA state={{ ...publicForm, ...main, id, source }} />
              </div>
              <div
                className={
                  displayTab === 'mutationalPattern' ? 'd-block' : 'd-none'
                }
              >
                <MutationalPattern
                  state={{ ...publicForm, ...main, id, source }}
                />
              </div>
              {source == 'user' && (
                <div
                  className={displayTab === 'cluster' ? 'd-block' : 'd-none'}
                >
                  <ClusteredIdentification
                    state={{ ...publicForm, ...main, params, id, source }}
                  />
                </div>
              )}
              <div className={displayTab === 'download' ? 'd-block' : 'd-none'}>
                <Download state={{ ...publicForm, ...main, id, source }} />
              </div>
            </>
          )}
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
