import { useCallback, useEffect, useMemo } from 'react';
import { Button, Nav, Form, Row, Col, Alert } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { useParams } from 'react-router-dom';
import { actions as extractionActions } from '@/services/store/extraction';
import { actions as modalActions } from '@/services/store/modal';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '@/components/controls/sidebar-container/sidebar-container';
import ExtractionForm from './extraction-form';
import Instructions from './instructions';
import Status from './status';
import SignatureMap from './signatureMap';
import TMB from '@/components/pages/exploration/tmb/tmb';
import TmbSignature from '@/components/pages/exploration/tmbSignature/tmbSignature';
import MsBurden from '@/components/pages/exploration/msBurden/msBurden';
import MsDecomposition from '@/components/pages/exploration/msDecomposition/msDecomposition';
import MsAssociation from '@/components/pages/exploration/msAssociation/msAssociation';
import MsLandscape from '@/components/pages/exploration/msLandscape/msLandscape';
import MsPrevalence from '@/components/pages/exploration/msPrevalence/msPrevalence';
import MsIndividual from '@/components/pages/exploration/msIndividual/msIndividual';
import { useRefreshQuery } from './apiSlice';
import { LoadingOverlay } from '@/components/controls/loading-overlay/loading-overlay';
import TabPanel from '@/components/controls/tab/tab';

const actions = { ...extractionActions, ...modalActions };

export default function Extraction() {
  const env = import.meta.env;
  const formLimits = env?.EXTRACTION_FORM_LIMIT
    ? JSON.parse(env.EXTRACTION_FORM_LIMIT)
    : {
        minimum_signatures: [1, 15],
        maximum_signatures: [1, 15],
        nmf_replicates: [1, 100],
        min_nmf_iterations: [1, 10000],
        max_nmf_iterations: [1, 1000000],
        nmf_test_conv: [1, 10000],
      };

  const dispatch = useDispatch();
  const mergeState = (state) => dispatch(actions.mergeExtraction(state));
  const { displayTab, openSidebar, explorationType, ...state } = useSelector(
    (state) => state.extraction
  );
  const id = useParams().id || state.id || false;

  const {
    data: refreshStatus,
    error,
    refetch: refreshExtraction,
  } = useRefreshQuery(id, { skip: !id });

  const status = refreshStatus?.status;
  const manifest = refreshStatus?.manifest;
  const params = refreshStatus?.params;

  const isDone = ['COMPLETED', 'FAILED'].includes(status?.status);
  const explorationId = (() => {
    if (isDone) {
      switch (explorationType) {
        case 'denovo':
          return manifest?.denovoId;
        case 'decomposed':
          return manifest?.decomposedId;
        default:
          return false;
      }
    } else {
      return false;
    }
  })();
  const id2 = explorationId
    ? explorationId === manifest.denovoId
      ? manifest.decomposedId
      : manifest.denovoId
    : false;

  const signatureMapState = useMemo(
    () => ({ id, params, manifest }),
    [id, params, manifest]
  );
  const explorationState = useMemo(
    () => ({ id: explorationId }),
    [explorationId]
  );
  const associationState = useMemo(
    () => ({ id: explorationId, id2 }),
    [explorationId, id2]
  );

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
    if (isDone || error) clearInterval(interval);
    return () => clearInterval(interval);
  }, [isDone, error, refreshState]);

  useEffect(() => {
    if (status && status.status === 'COMPLETED')
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
      disabled:
        explorationId === false ||
        !['default', 'SBS96', 'DBS72', 'ID83', 'RS32', 'SBS192'].includes(
          params?.args.context_type
        ),
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
          <ExtractionForm formLimits={formLimits} />
        </SidebarPanel>
        <MainPanel>
          {error && <Alert variant="danger">Results expired</Alert>}
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
                        disabled={displayTab === 'msAssociation'}
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
                    {/* <p>{manifest?.matrixFile}</p> */}
                    <p>
                      {manifest.matrixFile
                        ? manifest.matrixFile
                        : id.toLowerCase().includes('example')
                        ? 'example-data'
                        : 'public-data'}
                    </p>
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
                An error occurred during calculation:
                {status?.error?.message || 'INTERNAL ERROR'}. Please contact the
                site administrator for assistance if this issue persists.
              </p>
            </div>
          )}

          <TabPanel active={displayTab === 'instructions'}>
            <Instructions formLimits={formLimits} />
          </TabPanel>
          <TabPanel active={displayTab === 'status'}>
            <Status setDisplayTab={(tab) => mergeState({ displayTab: tab })} />
          </TabPanel>
          {status && status.status === 'COMPLETED' && (
            <>
              <TabPanel active={displayTab === 'signatureMap'}>
                <SignatureMap state={signatureMapState} />
              </TabPanel>
              <TabPanel active={displayTab === 'tmb'}>
                <TMB state={explorationState} />
              </TabPanel>
              <TabPanel active={displayTab === 'tmbSig'}>
                <TmbSignature state={explorationState} />
              </TabPanel>
              <TabPanel active={displayTab === 'msBurden'}>
                <MsBurden state={explorationState} />
              </TabPanel>
              <TabPanel active={displayTab === 'msDecomposition'}>
                <MsDecomposition state={explorationState} />
              </TabPanel>
              <TabPanel active={displayTab === 'msAssociation'}>
                <MsAssociation state={associationState} />
              </TabPanel>
              <TabPanel active={displayTab === 'msLandscape'}>
                <MsLandscape state={explorationState} />
              </TabPanel>
              <TabPanel active={displayTab === 'msPrevalence'}>
                <MsPrevalence state={explorationState} />
              </TabPanel>
              <TabPanel active={displayTab === 'msIndividual'}>
                <MsIndividual state={explorationState} />
              </TabPanel>
            </>
          )}
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
