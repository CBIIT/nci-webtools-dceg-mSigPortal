import React, { useEffect } from 'react';
import { Form, Row, Col } from 'react-bootstrap';
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
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import './visualization.scss';

const { Group, Label, Check } = Form;

export default function Visualize({ match }) {
  const { openSidebar, loading, source, submitted } = useSelector(
    (state) => state.visualize
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
        await fetch(`api/fetchResults/${id}`)
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
        await fetch(`api/fetchExample/${id}`)
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

  return (
    <div className="position-relative">
      <SidebarContainer
        className="m-3"
        collapsed={!openSidebar}
        onCollapsed={(e) => setOpenSidebar(!e)}
      >
        <SidebarPanel>
          <div className="p-3 bg-white border rounded">
            <Row>
              <Col sm="auto">
                <Group className="d-flex">
                  <Label className="mr-auto">
                    <h3 className="mb-2">Data Source</h3>
                  </Label>
                  <Check inline id="radioPublic" className="ml-4">
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
              <Col sm="auto" className="w-100">
                <UploadForm />
              </Col>
            </Row>
            <Row style={{ display: source == 'public' ? 'block' : 'none' }}>
              <Col sm="auto" className="w-100">
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
