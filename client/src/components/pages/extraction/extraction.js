import { useCallback, useEffect, useState } from 'react';
import { Button, Nav, Form } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { useParams } from 'react-router-dom';
import { actions as extractionActions } from '../../../services/store/extraction';
import { actions as modalActions } from '../../../services/store/modal';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '../../controls/sidebar-container/sidebar-container';
import ExtractionForm from './extraction-form';
import Instructions from './instructions';
import TMB from '../exposure/tmb/tmb';
import TmbSignature from '../exposure/tmbSignature/tmbSignature';
import MsBurden from '../exposure/msBurden/msBurden';
import MsDecomposition from '../exposure/msDecomposition/msDecomposition';
import MsAssociation from '../exposure/msAssociation/msAssociation';
import MsLandscape from '../exposure/msLandscape/msLandscape';
import MsPrevalence from '../exposure/msPrevalence/msPrevalence';
import { useStatusQuery, useManifestQuery } from './apiSlice';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

const actions = { ...extractionActions, ...modalActions };

export default function Extraction() {
  const dispatch = useDispatch();
  const mergeState = (state) => dispatch(actions.mergeExtraction(state));
  const [explorationType, setExplorationType] = useState('denovo');
  const { displayTab, openSidebar, ...state } = useSelector(
    (state) => state.extraction
  );
  const id = useParams().id || state.id || false;
  const { data: status, refetch: refetchStatus } = useStatusQuery(id, {
    skip: !id,
  });
  const { data: manifest, refetch: refetchManifest } = useManifestQuery(id, {
    skip: !id,
  });

  const isDone = ['COMPLETED', 'FAILED'].includes(status?.status);
  const explorationId = (() => {
    if (isDone) {
      if (explorationType === 'denovo') return manifest?.denovoId;
      else if (explorationType === 'decomposed') return manifest?.decomposedId;
      else return false;
    } else return false;
  })();

  const refreshState = useCallback(() => {
    refetchStatus();

    refetchManifest();
  }, [refetchStatus, refetchManifest]);

  useEffect(() => {
    const interval = setInterval(refreshState, 1000 * 60);
    if (isDone) clearInterval(interval);
    return () => clearInterval(interval);
  }, [isDone, refreshState]);

  useEffect(() => {
    if (status && status.status === 'COMPLETED')
      mergeState({ displayTab: 'tmb' });
  }, [isDone]);

  const tabs = [
    {
      id: 'instructions',
      name: 'Instructions',
    },
    {
      id: 'tmb',
      name: 'TMB',
    },
    {
      id: 'tmbSig',
      name: 'TMB Signatures',
    },
    {
      id: 'msBurden',
      name: 'MS Burden',
    },
    {
      id: 'msDecomposition',
      name: 'MS Decomposition',
    },
    {
      id: 'msAssociation',
      name: 'MS Association',
    },
    {
      id: 'msLandscape',
      name: 'MS Landscape',
    },
    {
      id: 'msPrevalence',
      name: 'MS Prevalence',
    },
    {
      id: 'msIndividual',
      name: 'MS Individual',
    },
  ];

  function handleExplorationType(e) {
    setExplorationType(e.target.value);
  }

  return (
    <div className="position-relative">
      <div className="mx-3">
        <div className="mx-3 bg-white border border-top-0">
          {/* for desktops and tablets */}
          <div className="d-none d-md-block">
            <Nav defaultActiveKey="instructions">
              {tabs.map(({ name, id }, i) => (
                <div key={name} className="d-inline-block">
                  <Button
                    variant="link"
                    className={`secondary-navlinks px-3 py-1 d-inline-block border-0 text-exploration rounded-0 ${
                      id === displayTab
                        ? 'bg-exploration text-white'
                        : 'text-exploration'
                    }`}
                    disabled={!explorationId}
                    style={{
                      textDecoration: 'none',
                      fontSize: '12pt',
                      color: '#837244',
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
          <div className="row d-md-none">
            <Nav defaultActiveKey="instructions">
              {tabs.map(({ name, id }, i) => {
                if (name)
                  return (
                    <div key={name} className="col-12 text-center">
                      <Button
                        variant="link"
                        className={
                          id === displayTab
                            ? 'secondary-navlinks px-3 py-1 d-inline-block border-0 bg-exploration text-white rounded-0'
                            : 'secondary-navlinks px-3 py-1 d-inline-block border-0 rounded-0'
                        }
                        style={{
                          textDecoration: 'none',
                          fontSize: '12pt',
                          color: '#837244',
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
        onCollapsed={(e) => mergeState({ openSidebar: !e })}
      >
        <SidebarPanel>
          {status && status.status === 'COMPLETED' && (
            <div className="p-3 bg-white border rounded mb-3">
              <Form.Group controlId="explorationType">
                <Form.Label>Exploration Calculation</Form.Label>
                <Form.Control as="select" onChange={handleExplorationType}>
                  <option value="denovo">Denovo</option>
                  <option value="decomposed">Decomposed</option>
                </Form.Control>
              </Form.Group>

              <Form.Group controlId="explorationExposure">
                <Form.Label>Exposure File</Form.Label>
                <Form.File
                  value={''} // set dummy value for file input
                  disabled={true}
                  label={
                    explorationType === 'denovo'
                      ? manifest?.denovoExposureFile
                      : manifest?.decomposedExposureFile
                  }
                  feedback="Please upload a data file"
                  custom
                />
              </Form.Group>
              <Form.Group controlId="explorationMatrix">
                <Form.Label>Matrix File</Form.Label>
                <Form.File
                  value={''} // set dummy value for file input
                  disabled={true}
                  label={manifest?.matrixFile}
                  feedback="Please upload a data file"
                  custom
                />
              </Form.Group>
              <Form.Group controlId="explorationSignature">
                <Form.Label>Signature File</Form.Label>
                <Form.File
                  value={''} // set dummy value for file input
                  disabled={true}
                  label={
                    explorationType === 'denovo'
                      ? manifest?.denovoSignatureFile
                      : manifest?.decomposedSignatureFile
                  }
                  feedback="Please upload a data file"
                  custom
                />
              </Form.Group>
            </div>
          )}
          <ExtractionForm setExplorationType={setExplorationType} />
        </SidebarPanel>
        <MainPanel>
          {status && status.status === 'IN_PROGRESS' && (
            <div className="border rounded bg-white mb-3 p-3">
              <p>
                Your analysis is currently in progress.
                {/* {params.sendNotification && (
                  <span> You will receive an email once it is complete.</span>
                )} */}
              </p>
              <LoadingOverlay active={true} />
            </div>
          )}
          {status && status.status === 'FAILED' && (
            <div className="border rounded bg-white mb-3 p-3">
              <p>
                Your analysis failed with the following error:{' '}
                {status?.error?.message || 'INTERNAL ERROR'}. Please contact the
                site administrator for assistance if this issue persists.
              </p>
            </div>
          )}
          <div className={displayTab === 'instructions' ? 'd-block' : 'd-none'}>
            <Instructions />
          </div>
          {status && status.status === 'COMPLETED' && (
            <>
              <div className={displayTab === 'tmb' ? 'd-block' : 'd-none'}>
                <TMB state={{ id: explorationId }} />
              </div>
              <div className={displayTab === 'tmbSig' ? 'd-block' : 'd-none'}>
                <TmbSignature state={{ id: explorationId }} />
              </div>
              <div className={displayTab === 'msBurden' ? 'd-block' : 'd-none'}>
                <MsBurden state={{ id: explorationId }} />
              </div>
              <div
                className={
                  displayTab === 'msDecomposition' ? 'd-block' : 'd-none'
                }
              >
                <MsDecomposition state={{ id: explorationId }} />
              </div>
              <div
                className={
                  displayTab === 'msAssociation' ? 'd-block' : 'd-none'
                }
              >
                <MsAssociation state={{ id: explorationId }} />
              </div>
              <div
                className={displayTab === 'msLandscape' ? 'd-block' : 'd-none'}
              >
                <MsLandscape state={{ id: explorationId }} />
              </div>
              <div
                className={displayTab === 'msPrevalence' ? 'd-block' : 'd-none'}
              >
                <MsPrevalence state={{ id: explorationId }} />
              </div>
              <div
                className={
                  displayTab === 'msIndivexplorationIdual'
                    ? 'd-block'
                    : 'd-none'
                }
              >
                {/* <MsIndividual state={{ id: explorationId }}  /> */}
                Under Construction
              </div>
            </>
          )}
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
