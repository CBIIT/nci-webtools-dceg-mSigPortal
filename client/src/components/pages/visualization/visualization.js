import React, { useEffect } from 'react';
import { Form, Row, Col, Nav, Button } from 'react-bootstrap';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '../../controls/sidebar-container/sidebar-container';
import UploadForm from './uploadForm';
import PublicForm from './publicForm';
import Results from './results';
import { useSelector } from 'react-redux';
import {
  dispatchVisualize,
  dispatchVisualizeResults,
  dispatchProfilerSummary,
  dispatchMutationalPattern,
  dispatchCosineSimilarity,
  dispatchProfileComparison,
  dispatchPCA,
  dispatchError,
  store,
  updateVisualize,
} from '../../../services/store';
import './visualization.scss';

const { Group, Label, Check } = Form;

export default function Visualize({ match }) {
  const { openSidebar, loading, source, submitted } = useSelector(
    (state) => state.visualize
  );
  const { displayTab, svgList } = useSelector(
    (state) => state.visualizeResults
  );
  const { type, id } = match.params;

  // when retrieving queued result, update id in store
  useEffect(() => {
    if (type == 'queue') {
      if (id) loadQueueResult(id);
    } else if (type == 'example') {
      if (id) loadExample(id);
    }
  }, [id]);

  function setOpenSidebar(bool) {
    dispatchVisualize({ openSidebar: bool });
  }

  async function loadQueueResult(id) {
    dispatchVisualize({
      loading: {
        active: true,
        content: 'Loading Queued Result',
        showIndicator: true,
      },
    });
    try {
      const { args, state, timestamp } = await (
        await fetch(`api/getQueueResults/${id}`)
      ).json();
      dispatchVisualize(state.visualize);
      dispatchVisualizeResults({ projectID: id });
    } catch (error) {
      dispatchError(error.toString());
    }
    dispatchVisualize({
      loading: { active: false },
    });
  }

  async function loadExample(id) {
    dispatchVisualize({
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
      dispatchVisualize(state.visualize);
      // rehydrate state if available
      if (state.profilerSummary) dispatchProfilerSummary(state.profilerSummary);
      if (state.mutationalPattern)
        dispatchMutationalPattern(state.mutationalPattern);
      if (state.cosineSimilarity)
        dispatchCosineSimilarity(state.cosineSimilarity);
      if (state.profileComparison)
        dispatchProfileComparison(state.profileComparison);
      if (state.pca) dispatchPCA(state.pca);
      if (state.visualizeResults) {
        dispatchVisualizeResults({
          ...state.visualizeResults,
          projectID: projectID,
        });
      } else dispatchVisualizeResults({ projectID: projectID });
    } catch (error) {
      dispatchError(error.toString());
    }
    dispatchVisualize({
      loading: { active: false },
    });
  }

  const subModules = [
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
              {subModules.map(({ name, id }) => (
                <div key={id} className="d-inline-block">
                  <Button
                    variant="link"
                    className={`secondary-navlinks px-3 py-1 d-inline-block border-0 ${
                      id == displayTab && svgList.length
                        ? 'active-secondary-navlinks'
                        : ''
                    }`}
                    active={id == displayTab && svgList.length}
                    disabled={!svgList.length}
                    style={{
                      textDecoration: 'none',
                      fontSize: '11pt',
                      color: 'black',
                      fontWeight: '500',
                    }}
                    onClick={() => dispatchVisualizeResults({ displayTab: id })}
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
              {subModules.map(({ name, id }) => (
                <div key={id} className="col-12 text-center">
                  <Button
                    variant="link"
                    className={
                      id == displayTab && svgList.length
                        ? 'secondary-navlinks px-3 py-1 d-inline-block border-0 active-secondary-navlinks'
                        : 'secondary-navlinks px-3 py-1 d-inline-block border-0'
                    }
                    style={{
                      textDecoration: 'none',
                      fontSize: '11pt',
                      color: 'black',
                      fontWeight: '500',
                    }}
                    onClick={() => dispatchVisualizeResults({ displayTab: id })}
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
                      onChange={(e) =>
                        store.dispatch(updateVisualize({ source: 'public' }))
                      }
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
                      onChange={(e) =>
                        store.dispatch(updateVisualize({ source: 'user' }))
                      }
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
                <UploadForm />
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
