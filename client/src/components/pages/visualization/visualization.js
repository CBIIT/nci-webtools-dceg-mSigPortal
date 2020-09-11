import React from 'react';
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
import { store, updateVisualize } from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import './visualization.scss';

const { Group, Label, Check } = Form;

export default function Visualize() {
  const { openSidebar, loading, source, submitted } = useSelector(
    (state) => state.visualize
  );

  function setOpenSidebar(bool) {
    store.dispatch(updateVisualize({ openSidebar: bool }));
  }

  return (
    <div className="position-relative">
      <SidebarContainer
        className="m-3"
        collapsed={!openSidebar}
        onCollapsed={(e) => setOpenSidebar(!e)}
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
                      onChange={(e) =>
                        store.dispatch(updateVisualize({ source: 'user' }))
                      }
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
                      onChange={(e) =>
                        store.dispatch(updateVisualize({ source: 'public' }))
                      }
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
          <div
            className="p-3 shadow-sm bg-white"
            style={{ minHeight: '420px' }}
          >
            <LoadingOverlay
              active={loading.active}
              content={loading.content}
              showIndicator={loading.showIndicator}
            />
            <Results setOpenSidebar={(e) => setOpenSidebar(e)} />
          </div>
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
