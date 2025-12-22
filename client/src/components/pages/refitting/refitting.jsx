import React, { useCallback, useEffect } from 'react';
import { Button, Nav, Alert } from 'react-bootstrap';
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
import Results from './results';
import RefittingForm from './refitting-form';
import { useRefittingStatusQuery } from './apiSlice';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

export default function Refitting() {
  const { displayTab, openSidebar, ...state } = useSelector(
    (state) => state.refitting.main
  );
  const dispatch = useDispatch();
  const urlParamId = useParams()?.id;
  const jobId = urlParamId || state.id || false;

  const {
    data: refreshStatus,
    error,
    refetch: refreshRefitting,
  } = useRefittingStatusQuery(jobId, { skip: !jobId });

  const status = refreshStatus?.status;
  const isDone = ['COMPLETED', 'FAILED'].includes(status?.status);

  const refreshState = useCallback(() => {
    refreshRefitting();
  }, [refreshRefitting]);

  useEffect(() => {
    const interval = setInterval(refreshState, 1000 * 60);
    if (isDone || error) clearInterval(interval);
    return () => clearInterval(interval);
  }, [isDone, error, refreshState]);

  // Save jobId from URL to store when job is valid
  useEffect(() => {
    if (urlParamId && status) {
      dispatch(
        refittingActions.mergeRefitting({
          main: { id: urlParamId, submitted: true },
        })
      );
    }
  }, [urlParamId, status, dispatch]);

  const setDisplayTab = (tab) =>
    dispatch(
      refittingActions.mergeRefitting({
        main: { displayTab: tab },
      })
    );

  const setOpenSidebar = (open) =>
    dispatch(
      refittingActions.mergeRefitting({
        main: { openSidebar: open },
      })
    );

  useEffect(() => {
    if (status && status.status === 'COMPLETED') {
      setDisplayTab('results');
    }
  }, [status]);

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
      id: 'results',
      name: 'Results',
      disabled: !jobId || status?.status !== 'COMPLETED',
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
          {error && <Alert variant="danger">Results expired</Alert>}
          {status?.status === 'SUBMITTED' && (
            <div className="border rounded bg-white mb-3 p-3">
              <p>
                Your job has been submitted. You will receive an email once it
                is complete.
              </p>
            </div>
          )}
          {status?.status === 'IN_PROGRESS' && (
            <div className="border rounded bg-white mb-3 p-3">
              <p>Your analysis is currently in progress.</p>
              <LoadingOverlay active={true} />
            </div>
          )}
          {status?.status === 'FAILED' && (
            <div className="border rounded bg-white mb-3 p-3">
              <p>
                An error occurred during calculation:{' '}
                {status?.error || 'INTERNAL ERROR'}. Please contact the site
                administrator for assistance if this issue persists.
              </p>
            </div>
          )}

          <div className={displayTab === 'instructions' ? 'd-block' : 'd-none'}>
            <Instructions />
          </div>
          <div className={displayTab === 'status' ? 'd-block' : 'd-none'}>
            <Status setDisplayTab={setDisplayTab} />
          </div>
          {status && status.status === 'COMPLETED' && (
            <>
              <div className={displayTab === 'results' ? 'd-block' : 'd-none'}>
                <Results jobId={jobId} />
              </div>
            </>
          )}
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
