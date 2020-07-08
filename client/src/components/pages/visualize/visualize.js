import React, { useState } from 'react';
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
    <div className="container">
      <SidebarContainer
        className="my-3"
        collapsed={!openSidebar}
        onCollapsed={(e) => setOpenSidebar(!e)}
      >
        <SidebarPanel className="col-lg-4">
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
        </SidebarPanel>
        <MainPanel className="col-lg-8">
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
