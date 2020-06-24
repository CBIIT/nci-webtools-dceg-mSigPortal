import React, { useState } from 'react';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '../../controls/sidebar-container/sidebar-container';
import UploadForm from './uploadForm';
import Results from './results';
import './visualize.scss';

export default function Visualize() {
  const [openSidebar, setOpenSidebar] = useState(true);

  return (
    <div className="container">
      <SidebarContainer
        className="my-3"
        collapsed={!openSidebar}
        onCollapsed={(collapsed) => setOpenSidebar(!collapsed)}
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
                <UploadForm />
              </div>
            </div>
          </div>
        </SidebarPanel>
        <MainPanel className="col-lg-8">
          <div
            className="p-3 shadow-sm bg-white"
            style={{ minHeight: '420px' }}
          >
            <Results />
          </div>
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
