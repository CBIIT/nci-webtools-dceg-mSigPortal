import React from 'react';
import { Form, Row, Col, Nav, Button } from 'react-bootstrap';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '../../controls/sidebar-container/sidebar-container';
import Univariate from './univariate';
import UserForm from './userForm';
import PublicForm from './publicForm';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../services/store/association';
import { actions as modalActions } from '../../../services/store/modal';
import './association.scss';

const actions = { ...visualizationActions, ...modalActions };
const { Group, Label, Check } = Form;

export default function Association() {
  const dispatch = useDispatch();
  const mergeState = async (state) =>
    await dispatch(actions.mergeAssociation({ associationState: state }));

  const {
    displayTab,
    openSidebar,
    submitted,
    expVarList,
    source,
    assocVarData,
    loadingData,
  } = useSelector((state) => state.association.associationState);

  const tabs = [
    { name: 'Univariate', id: 'univariate' },
    { name: 'Multivariate', id: 'multivariate' },
  ];

  return (
    <div className="position-relative">
      <div className="mx-3">
        <div className="mx-3 bg-white border border-top-0">
          {/* for desktops and tablets */}
          <div className="d-none d-md-block">
            <Nav defaultActiveKey="profilerSummary">
              {tabs.map(({ name, id }) => (
                <div key={id} className="d-inline-block">
                  <Button
                    variant="link"
                    className={`secondary-navlinks px-3 py-1 d-inline-block border-0 ${
                      id == displayTab && Object.keys(expVarList).length
                        ? 'active-secondary-navlinks'
                        : ''
                    }`}
                    active={id == displayTab && Object.keys(expVarList).length}
                    disabled={!Object.keys(expVarList).length}
                    style={{
                      textDecoration: 'none',
                      fontSize: '12pt',
                      color: 'black',
                      fontWeight: '500',
                    }}
                    onClick={() => mergeState({ displayTab: id })}
                  >
                    {name}
                  </Button>
                  <div className="d-md-none w-100"></div>
                </div>
              ))}
            </Nav>
          </div>
          {/* for mobile devices */}
          <div className="row d-md-none">
            <Nav defaultActiveKey="summary">
              {tabs.map(({ name, id }) => (
                <div key={id} className="col-12 text-center">
                  <Button
                    variant="link"
                    className={
                      id == displayTab && Object.keys(expVarList).length
                        ? 'secondary-navlinks px-3 py-1 d-inline-block border-0 active-secondary-navlinks'
                        : 'secondary-navlinks px-3 py-1 d-inline-block border-0'
                    }
                    style={{
                      textDecoration: 'none',
                      fontSize: '12pt',
                      color: 'black',
                      fontWeight: '500',
                    }}
                    onClick={() => mergeState({ displayTab: id })}
                  >
                    {name}
                  </Button>
                  <div className="d-md-none w-100"></div>
                </div>
              ))}
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
                      disabled={submitted}
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
                      disabled={submitted}
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
            <Row
              style={{
                display: source == 'user' ? 'block' : 'none',
              }}
            >
              <Col lg="auto" className="w-100">
                {/* <UserForm /> */}
              </Col>
            </Row>
            <Row style={{ display: source == 'public' ? 'block' : 'none' }}>
              <Col lg="auto" className="w-100">
                <PublicForm />
              </Col>
            </Row>
          </div>
          <hr className="d-lg-none" style={{ opacity: 0 }}></hr>
        </SidebarPanel>
        <MainPanel>
          <>
            {assocVarData.length ? (
              displayTab == 'univariate' ? (
                <Univariate />
              ) : (
                <></>
              )
            ) : (
              <div className="bg-white border rounded py-3 px-4">
                <LoadingOverlay active={loadingData} />
                <h4>Instructions</h4>
                <p>
                  Choose a Data Source and its associated options to submit a
                  query using the panel on the left
                </p>
                <hr />
                <h4>Data Source</h4>
                <p>
                  Public: Perform analysis using data available on the website
                </p>
                <p>User: Upload your own data</p>
              </div>
            )}
          </>
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
