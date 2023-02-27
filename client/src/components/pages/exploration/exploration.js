import React, { useState, useEffect } from 'react';
import { Form, Row, Col, Button, Nav } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { saveAs } from 'file-saver';
import { useParams } from 'react-router-dom';
import PublicForm from './publicForm/publicForm';
import UserForm from './userForm/userForm';
import Instructions from './instructions';
import TMB from './tmb/tmb.js';
import TmbSig from './tmbSignature/tmbSignature.js';
import MsBurden from './msBurden/msBurden.js';
import MsAssociation from './msAssociation/msAssociation.js';
import MsDecomposition from './msDecomposition/msDecomposition.js';
import MsLandscape from './msLandscape/msLandscape.js';
import MsPrevalence from './msPrevalence/msPrevalence.js';
import MsIndividual from './msIndividual/msIndividual.js';
import Download from './download';
import { actions as exposureActions } from '../../../services/store/exploration';
import { actions as modalActions } from '../../../services/store/modal';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '../../controls/sidebar-container/sidebar-container';

const actions = { ...exposureActions, ...modalActions };
const { Group, Label, Check } = Form;

export default function Exploration() {
  const dispatch = useDispatch();

  const mergeState = (state) =>
    dispatch(actions.mergeExploration({ main: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const { publicForm, main } = useSelector((state) => state.exploration);

  const { exampleName } = useParams();

  // const [variableFileObj, setVariable] = useState(new File([], ''));

  const { displayTab, source, loading, id, openSidebar, submitted } =
    main;

  // load example if available
  useEffect(() => {
    if (exampleName) loadExample(exampleName);
  }, [exampleName]);

  async function loadExample(id) {
    mergeState({
      loading: {
        active: true,
        // content: 'Loading Example',
        // showIndicator: true,
      },
    });
    try {
      const { state } = await (
        await fetch(`web/getExposureExample/${id}`)
      ).json();

      dispatch(actions.mergeExploration(state));
    } catch (error) {
      mergeError('Example does not exist');
    }
    mergeState({
      loading: false,
      submitted: true,
      openSidebar: false,
    });
  }

  async function submitR(fn, args, id = id) {
    try {
      const response = await fetch(`web/explorationWrapper`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          fn: fn,
          args: args,
          id,
        }),
      });

      if (response.ok) {
        return await response.json();
      } else {
        mergeError('R submit failed');
      }
    } catch (err) {
      mergeError(err.message);
    }
  }

  // when using public data and only need to upload a variable data file
  // also used to create a work directory and id
  // async function uploadVariable() {
  //   return new Promise(async (resolve, reject) => {
  //     try {
  //       const data = new FormData();
  //       if (variableFileObj.size) data.append('variableFile', variableFileObj);

  //       let response = await fetch(`web/upload`, {
  //         method: 'POST',
  //         body: data,
  //       });

  //       if (!response.ok) {
  //         const { msg, error } = await response.json();

  //         mergeError([msg, error]);
  //         reject(error);
  //       } else {
  //         const { id } = await response.json();
  //         await mergeState({ id });
  //         resolve(id);
  //       }
  //     } catch (err) {
  //       mergeError(err.message);
  //       reject(err);
  //     }
  //   });
  // }

  // function handleVariable(file) {
  //   console.log(file);
  //   setVariable(file);
  //   mergeMsLandscape({ variableFile: file.name });
  // }

  async function exposureDownload() {
    try {
      const { output, id, stdout } = await submitR('exposureDownload', {
        study: publicForm.study.value,
        strategy: publicForm.strategy.value,
        rsSet: publicForm.signatureSet.value,
        cancerType: publicForm.cancer.value,
      });

      const file = await fetch(`web/data/${output.path}`);
      if (file.ok) {
        saveAs(await file.blob(), output.filename);
      } else {
        mergeError(`public data is not available`);
      }
    } catch (err) {
      console.log(err);
    }
  }

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
    source == 'public' ? (
      {
        id: 'download',
        name: 'Download',
      }
    ) : (
      <></>
    ),
  ];

  return (
    <div className="position-relative">
      <div className="mx-3">
        <div className="mx-3 bg-white border border-top-0">
          {/* for desktops and tablets */}
          <div className="d-none d-md-block">
            <Nav defaultActiveKey="profilerSummary">
              {tabs.map(({ name, id }) => {
                if (name)
                  return (
                    <div key={id} className="d-inline-block">
                      <Button
                        variant="link"
                        className={`secondary-navlinks px-3 py-1 d-inline-block border-0 text-exploration rounded-0 ${
                          id == displayTab
                            ? 'bg-exploration text-white'
                            : 'text-exploration'
                        }`}
                        active={id == displayTab && submitted}
                        disabled={
                          id != 'instructions' &&
                          !(source == 'public' ? submitted : main.id)
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
                    </div>
                  );
              })}
            </Nav>
          </div>

          {/* for mobile devices */}
          <div className="row d-md-none">
            <Nav defaultActiveKey="summary">
              {tabs.map(({ name, id }) => {
                if (name)
                  return (
                    <div key={id} className="col-12 text-center">
                      <Button
                        variant="link"
                        className={
                          id == displayTab
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
          <div className="p-3 bg-white border rounded">
            <Row>
              <Col lg="auto">
                <Group>
                  <Label className="mr-4">Data Source</Label>
                  <Check inline id="radioPublic">
                    <Check.Input
                      disabled={loading || submitted || submitted}
                      type="radio"
                      value="public"
                      checked={source == 'public'}
                      onChange={(e) => mergeState({ source: 'public' })}
                    />
                    <Check.Label className="font-weight-normal">
                      Public
                    </Check.Label>
                  </Check>
                  <Check inline id="radioUser">
                    <Check.Input
                      disabled={loading || submitted || submitted}
                      type="radio"
                      value="user"
                      checked={source == 'user'}
                      onChange={(e) => mergeState({ source: 'user' })}
                    />
                    <Check.Label className="font-weight-normal">
                      User
                    </Check.Label>
                  </Check>
                </Group>
              </Col>
            </Row>
            <Row>
              <Col lg="12" className="w-100">
                {source == 'public' ? <PublicForm /> : <UserForm />}
              </Col>
            </Row>
          </div>
        </SidebarPanel>
        <MainPanel>
          <div className={displayTab == 'instructions' ? 'd-block' : 'd-none'}>
            <Instructions />
          </div>
          <div className={displayTab == 'tmb' ? 'd-block' : 'd-none'}>
            <TMB state={{ ...publicForm, ...main }} />
          </div>
          <div className={displayTab == 'tmbSig' ? 'd-block' : 'd-none'}>
            <TmbSig state={{ ...publicForm, ...main }} />
          </div>
          <div className={displayTab == 'msBurden' ? 'd-block' : 'd-none'}>
            <MsBurden state={{ ...publicForm, ...main }} />
          </div>
          <div
            className={displayTab == 'msDecomposition' ? 'd-block' : 'd-none'}
          >
            <MsDecomposition state={{ ...publicForm, ...main }} />
          </div>
          <div className={displayTab == 'msAssociation' ? 'd-block' : 'd-none'}>
            <MsAssociation state={{ ...publicForm, ...main }} />
          </div>
          <div className={displayTab == 'msLandscape' ? 'd-block' : 'd-none'}>
            <MsLandscape state={{ ...publicForm, ...main }} />
          </div>
          <div className={displayTab == 'msPrevalence' ? 'd-block' : 'd-none'}>
            <MsPrevalence state={{ ...publicForm, ...main }} />
          </div>
          <div className={displayTab == 'msIndividual' ? 'd-block' : 'd-none'}>
            <MsIndividual state={{ ...publicForm, ...main }} />
          </div>
          <div
            className={
              displayTab == 'download' && source == 'public'
                ? 'd-block'
                : 'd-none'
            }
          >
            <Download exposureDownload={exposureDownload} />
          </div>
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
