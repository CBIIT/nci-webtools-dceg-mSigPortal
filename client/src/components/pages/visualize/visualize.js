import React from 'react';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '../../controls/sidebar-container/sidebar-container';
import UploadForm from './uploadForm';
import Results from './results';
import { useSelector } from 'react-redux';
import { store, updateVisualize } from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import './visualize.scss';

export default function Visualize() {
  const { openSidebar, loading } = useSelector((state) => state.visualize);

  function setOpenSidebar(bool) {
    store.dispatch(updateVisualize({ openSidebar: bool }));
  }

  return (
    <div className="position-relative mx-4">
      <SidebarContainer
        className="m-3"
        collapsed={!openSidebar}
        onCollapsed={(e) => setOpenSidebar(!e)}
      >
        <SidebarPanel>
          <div className="p-3 shadow-sm bg-white">
            <div className="row">
              <div className="col-sm-auto">
                <h3 className="mb-2">Parameters</h3>
              </div>
            </div>
            <div className="row">
              <div className="col-sm-auto w-100">
                <UploadForm setOpenSidebar={(e) => setOpenSidebar(e)} />
              </div>
            </div>
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
            <Results />
          </div>
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
