import React, { useState, useEffect } from 'react';
import { Button, Nav } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { useParams } from 'react-router-dom';
import { actions as refittingActions } from '../../../services/store/refitting';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '../../controls/sidebar-container/sidebar-container';
import Instructions from './instructions';
import Status from './status';
import TargetedSequencing from './targetedSequencing';
import RefittingForm from './refitting-form';

export default function Refitting() {
  const { displayTab, openSidebar } = useSelector((state) => state.refitting.main);
  const dispatch = useDispatch();
  const { id: jobId } = useParams();
  
  const setDisplayTab = (tab) => dispatch(refittingActions.mergeRefitting({ 
    main: { displayTab: tab } 
  }));
  
  const setOpenSidebar = (open) => dispatch(refittingActions.mergeRefitting({ 
    main: { openSidebar: open } 
  }));

  // Auto-switch to targeted sequencing tab when job ID is provided
  useEffect(() => {
    if (jobId) {
      setDisplayTab('targetedSequencing');
    }
  }, [jobId]);

  const tabs = [
    {
      id: 'instructions',
      name: 'Instructions',
      disabled: false,
    },
    {
      id: 'status',
      name: 'Status',
      disabled: false,
    },
    {
      id: 'targetedSequencing',
      name: 'Targeted Sequencing',
      disabled: !jobId, // Only enabled when a job ID is present
    },
  ];

  return (
    <div className="position-relative refitting-page">
      <h1 className="sr-only">Refitting</h1>
      <div className="mx-3">
        <div className="mx-3 bg-white border border-top-0">
          {/* for desktops and tablets */}
          <div className="d-none d-md-block">
            <Nav defaultActiveKey="instructions">
              {tabs.map(({ name, id, disabled }, i) => (
                <div key={name} className="d-inline-block">
                  <Button
                    variant="link"
                    className={`secondary-navlinks px-3 py-1 d-inline-block border-0 text-refitting rounded-0 ${
                      id === displayTab
                        ? 'bg-refitting text-white'
                        : 'text-refitting'
                    }`}
                    disabled={disabled}
                    style={{
                      textDecoration: 'none',
                      fontSize: '12pt',
                      color: '#406024',
                      fontWeight: '500',
                    }}
                    onClick={() => setDisplayTab(id)}
                  >
                    {name}
                  </Button>
                </div>
              ))}
            </Nav>
          </div>

          {/* for mobile devices */}
          <div className="row d-md-none">
            <Nav defaultActiveKey="instructions">
              {tabs.map(({ name, id, disabled }, i) => {
                if (name)
                  return (
                    <div key={name} className="col-12 text-center">
                      <Button
                        variant="link"
                        className={
                          id === displayTab
                            ? 'secondary-navlinks px-3 py-1 d-inline-block border-0 bg-refitting text-white rounded-0'
                            : 'secondary-navlinks px-3 py-1 d-inline-block border-0 rounded-0'
                        }
                        disabled={disabled}
                        style={{
                          textDecoration: 'none',
                          fontSize: '12pt',
                          color: '#689f39',
                          fontWeight: '500',
                        }}
                        onClick={() => setDisplayTab(id)}
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
        onCollapsed={(e) => setOpenSidebar(!e)}
      >
        <SidebarPanel>
          <RefittingForm />
        </SidebarPanel>
        <MainPanel>
          <div className={displayTab === 'instructions' ? 'd-block' : 'd-none'}>
            <Instructions />
          </div>
          <div className={displayTab === 'status' ? 'd-block' : 'd-none'}>
            <Status />
          </div>
          <div className={displayTab === 'targetedSequencing' ? 'd-block' : 'd-none'}>
            <TargetedSequencing jobId={jobId} />
          </div>
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
