import React, { useEffect } from 'react';
import { Form, Row, Col, Nav, Button } from 'react-bootstrap';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '../../controls/sidebar-container/sidebar-container';
import UserForm from './userForm';
import PublicForm from './publicForm';
import Results from './results';
import './visualization.scss';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../services/store/visualization';
import { actions as modalActions } from '../../../services/store/modal';

const actions = { ...visualizationActions, ...modalActions };
const { Group, Label, Check } = Form;

export default function Visualize({ match }) {
  const dispatch = useDispatch();
  const visualization = useSelector((state) => state.visualization);

  const mergeVisualize = (state) =>
    dispatch(actions.mergeVisualization({ visualize: state }));
  const mergeResults = (state) =>
    dispatch(actions.mergeVisualization({ results: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const { openSidebar, loading, source, submitted } = visualization.visualize;
  const { displayTab, svgList } = visualization.results;
  const { type, id } = match.params;

  // when retrieving queued result, update id in store
  useEffect(() => {
    if (id && !loading.active && !submitted) {
      if (type == 'queue') {
        loadQueueResult(id);
      } else if (type == 'example') {
        loadExample(id);
      }
    }
  }, [id, loading.active]);

  function setOpenSidebar(bool) {
    mergeVisualize({ openSidebar: bool });
  }

  async function loadQueueResult(id) {
    mergeVisualize({
      loading: {
        active: true,
        content: 'Loading Queued Result',
        showIndicator: true,
      },
    });
    try {
      const { args, visualization, timestamp } = await (
        await fetch(`api/getQueueResults/${id}`)
      ).json();
      dispatch(actions.mergeVisualization(visualization));
    } catch (error) {
      mergeError(error.toString());
    }
    mergeVisualize({
      loading: { active: false },
    });
  }

  async function loadExample(id) {
    mergeVisualize({
      loading: {
        active: true,
        content: 'Loading Example',
        showIndicator: true,
      },
    });
    try {
      const { projectID, state } = await (
        await fetch(`api/getVisExample/${id}`)
      ).json();
      dispatch(
        actions.mergeVisualization({
          ...state,
          results: { ...state.results, projectID: projectID },
        })
      );
    } catch (error) {
      mergeError(error.toString());
    }
    mergeVisualize({
      loading: { active: false },
      submitted: true,
    });
  }

  const tabs = [
    { name: 'Profiler Summary', id: 'profilerSummary' },
    { name: 'Mutational Profiles', id: 'mutationalProfiles' },
    { name: 'Cosine Similarity', id: 'cosineSimilarity' },
    {
      name: 'Mutational Pattern Enrichment Analysis',
      id: 'mutationalPattern',
    },
    { name: 'Profile Comparison', id: 'profileComparison' },
    { name: 'PCA', id: 'pca' },
    { name: 'Kataegis Identification', id: 'kataegis' },
    { name: 'Download', id: 'download' },
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
                      id == displayTab && Object.keys(svgList).length
                        ? 'active-secondary-navlinks'
                        : ''
                    }`}
                    active={id == displayTab && Object.keys(svgList).length}
                    disabled={!Object.keys(svgList).length}
                    style={{
                      textDecoration: 'none',
                      fontSize: '11pt',
                      color: 'black',
                      fontWeight: '500',
                    }}
                    onClick={() => mergeResults({ displayTab: id })}
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
                      fontSize: '11pt',
                      color: 'black',
                      fontWeight: '500',
                    }}
                    onClick={() => mergeResults({ displayTab: id })}
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
        onCollapsed={(e) => setOpenSidebar(!e)}
      >
        <SidebarPanel>
          <div className="p-3 bg-white border rounded">
            <Row>
              <Col lg="auto">
                <Group>
                  <Label className="mr-4">Data Source</Label>
                  <Check inline id="radioPublic">
                    <Check.Input
                      disabled={submitted}
                      type="radio"
                      value="public"
                      checked={source == 'public'}
                      onChange={(e) => mergeVisualize({ source: 'public' })}
                    />
                    <Check.Label className="font-weight-normal">
                      Public
                    </Check.Label>
                  </Check>
                  <Check inline id="radioUser">
                    <Check.Input
                      disabled={submitted}
                      type="radio"
                      value="user"
                      checked={source == 'user'}
                      onChange={(e) => mergeVisualize({ source: 'user' })}
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
        <MainPanel>
          <Results setOpenSidebar={(e) => setOpenSidebar(e)} />
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
