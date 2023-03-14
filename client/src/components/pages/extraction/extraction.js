import { useCallback, useEffect } from 'react';
import { Button, Nav, Form, Row, Col } from 'react-bootstrap';
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
import Status from './status';
import SignatureMap from './signatureMap';
import TMB from '../exploration/tmb/tmb';
import TmbSignature from '../exploration/tmbSignature/tmbSignature';
import MsBurden from '../exploration/msBurden/msBurden';
import MsDecomposition from '../exploration/msDecomposition/msDecomposition';
import MsAssociation from '../exploration/msAssociation/msAssociation';
import MsLandscape from '../exploration/msLandscape/msLandscape';
import MsPrevalence from '../exploration/msPrevalence/msPrevalence';
import MsIndividual from '../exploration/msIndividual/msIndividual';
import { useRefreshQuery } from './apiSlice';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

const actions = { ...extractionActions, ...modalActions };

export default function Extraction() {
  const dispatch = useDispatch();
  const mergeState = (state) => dispatch(actions.mergeExtraction(state));
  const { displayTab, openSidebar, explorationType, ...state } = useSelector(
    (state) => state.extraction
  );
  const id = useParams().id || state.id || false;

  const { data: refreshStatus, refetch: refreshExtraction } = useRefreshQuery(
    id,
    { skip: !id }
  );

  const status = refreshStatus?.status;
  const manifest = refreshStatus?.manifest;
  const params = refreshStatus?.params;
  // const { data: status, refetch: refetchStatus } = useStatusQuery(id, {
  //   skip: !id,
  // });
  // const { data: manifest, refetch: refetchManifest } = useManifestQuery(id, {
  //   skip: !id,
  // });

  const isDone = ['COMPLETED', 'FAILED'].includes(status?.status);
  const explorationId = (() => {
    if (isDone) {
      if (explorationType === 'denovo') return manifest?.denovoId;
      else if (explorationType === 'decomposed') return manifest?.decomposedId;
      else return false;
    } else return false;
  })();

  const refreshState = useCallback(() => {
    refreshExtraction();
  }, [refreshExtraction]);
  // const refreshState = useCallback(() => {
  //   refreshExtraction();
  //   refetchStatus();
  //   refetchManifest();
  // }, [refreshExtraction, refetchStatus, refetchManifest]);

  useEffect(() => {
    const interval = setInterval(refreshState, 1000 * 60);
    if (isDone) clearInterval(interval);
    return () => clearInterval(interval);
  }, [isDone, refreshState]);

  useEffect(() => {
    if (status && status.status === 'COMPLETED' && displayTab == 'instructions')
      mergeState({ displayTab: 'signatureMap', openSidebar: false });
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
      id: 'signatureMap',
      name: 'Signature Map',
      disabled: status?.status !== 'COMPLETED',
    },
    {
      id: 'tmb',
      name: 'TMB',
      disabled: explorationId === false,
    },
    {
      id: 'tmbSig',
      name: 'TMB Signatures',
      disabled: explorationId === false,
    },
    {
      id: 'msBurden',
      name: 'MS Burden',
      disabled: explorationId === false,
    },
    {
      id: 'msDecomposition',
      name: 'MS Decomposition',
      disabled: explorationId === false,
    },
    {
      id: 'msAssociation',
      name: 'MS Association',
      disabled: explorationId === false,
    },
    {
      id: 'msLandscape',
      name: 'MS Landscape',
      disabled: explorationId === false,
    },
    {
      id: 'msPrevalence',
      name: 'MS Prevalence',
      disabled: explorationId === false,
    },
    {
      id: 'msIndividual',
      name: 'MS Individual',
      disabled: explorationId === false,
    },
  ];

  function handleExplorationType(e) {
    mergeState({ explorationType: e.target.value });
  }

  return (
    <div className="position-relative">
      <div className="mx-3">
        <div className="mx-3 bg-white border border-top-0">
          {/* for desktops and tablets */}
          <div className="d-none d-md-block">
            <Nav defaultActiveKey="instructions">
              {tabs.map(({ name, id, disabled }, i) => (
                <div key={name} className="d-inline-block">
                  <Button
                    variant="link"
                    className={`secondary-navlinks px-3 py-1 d-inline-block border-0 text-exploration rounded-0 ${
                      id === displayTab
                        ? 'bg-extraction text-white'
                        : 'text-extraction'
                    }`}
                    disabled={disabled}
                    style={{
                      textDecoration: 'none',
                      fontSize: '12pt',
                      color: '#42688b',
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
              {tabs.map(({ name, id, disabled }, i) => {
                if (name)
                  return (
                    <div key={name} className="col-12 text-center">
                      <Button
                        variant="link"
                        className={
                          id === displayTab
                            ? 'secondary-navlinks px-3 py-1 d-inline-block border-0 bg-extraction text-white rounded-0'
                            : 'secondary-navlinks px-3 py-1 d-inline-block border-0 rounded-0'
                        }
                        disabled={disabled}
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
          <ExtractionForm />
        </SidebarPanel>
        <MainPanel>
          {status &&
            status.status === 'COMPLETED' &&
            !['instructions', 'status', 'signatureMap'].includes(
              displayTab
            ) && (
              <div className="p-3 bg-white border rounded mb-3">
                <Row>
                  <Col md="auto">
                    <Form.Group controlId="explorationType">
                      <Form.Label>Exploration Calculation</Form.Label>
                      <Form.Control
                        as="select"
                        value={explorationType}
                        onChange={handleExplorationType}
                      >
                        <option value="denovo">Denovo</option>
                        <option value="decomposed">Decomposed</option>
                      </Form.Control>
                    </Form.Group>
                  </Col>
                  <Col md="auto">
                    <b>Exposure File</b>
                    <p>
                      {explorationType === 'denovo'
                        ? manifest?.denovoExposureInput
                        : manifest?.decomposedExposureInput}
                    </p>
                  </Col>
                  <Col md="auto">
                    <b>Matrix File</b>
                    <p>{manifest?.matrixFile}</p>
                  </Col>
                  <Col md="auto">
                    <b>Signature File</b>
                    <p>
                      {explorationType === 'denovo'
                        ? manifest?.denovoSignatureInput
                        : manifest?.decomposedSignatureInput}
                    </p>
                  </Col>
                </Row>
              </div>
            )}

          {status && status.status === 'SUBMITTED' && (
            <div className="border rounded bg-white mb-3 p-3">
              <p>
                Your job has been submitted. You will receive an email once it
                is complete.
              </p>
            </div>
          )}
          {status && status.status === 'IN_PROGRESS' && (
            <div className="border rounded bg-white mb-3 p-3">
              <p>Your analysis is currently in progress.</p>
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
          <div className={displayTab === 'status' ? 'd-block' : 'd-none'}>
            <Status />
          </div>
          {status && status.status === 'COMPLETED' && (
            <>
              <div
                className={displayTab === 'signatureMap' ? 'd-block' : 'd-none'}
              >
                <SignatureMap state={{ id, params, manifest }} />
              </div>
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
                className={displayTab === 'msIndividual' ? 'd-block' : 'd-none'}
              >
                <MsIndividual state={{ id: explorationId }} />
              </div>
            </>
          )}
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
