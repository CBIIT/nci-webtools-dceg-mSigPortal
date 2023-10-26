import React from 'react';
import { Form, Row, Col, Nav, Button } from 'react-bootstrap';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '../../controls/sidebar-container/sidebar-container';
import Instructions from '../association/instructions';
import Univariable from './univariable';
import Multivariable from './multivariable';
// import UserForm from './userForm';
import PublicForm from './publicForm/publicForm';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../services/store/association';
import { actions as modalActions } from '../../../services/store/modal';
import './association.scss';

const actions = { ...visualizationActions, ...modalActions };
const { Group, Label, Check } = Form;

export default function Association() {
  const dispatch = useDispatch();
  const mergeState = async (state) =>
    await dispatch(actions.mergeAssociation({ main: state }));

  const {
    displayTab,
    openSidebar,
    submitted,
    expVarList,
    source,
    loadingData,
  } = useSelector((state) => state.association.main);

  const tabs = [
    {
      name: 'Instructions',
      id: 'instructions',
      component: <Instructions loading={loadingData} />,
    },
    { name: 'Univariable', id: 'univariable', component: <Univariable /> },
    {
      name: 'Multivariable',
      id: 'multivariable',
      component: <Multivariable />,
    },
  ];

  return (
    <div className="position-relative">
      <div className="mx-3">
        <div className="mx-3 bg-white border border-top-0">
          {/* for desktops and tablets */}
          <div className="d-none d-md-block">
            <Nav defaultActiveKey="profilerSummary">
              {tabs.map(({ name, id }) => (
                <div key={id} className="d-inline-block ">
                  <Button
                    variant="link"
                    className={`secondary-navlinks px-3 py-1 d-inline-block border-0 association rounded-0 ${
                      id == displayTab ? 'bg-association text-white' : ''
                    }`}
                    active={id == displayTab && submitted}
                    disabled={
                      id != 'instructions' && !Object.keys(expVarList).length
                    }
                    style={{
                      textDecoration: 'none',
                      fontSize: '12pt',
                      color: '#b83d47',
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
            <Nav defaultActiveKey="summary">
              {tabs.map(({ name, id }) => (
                <div key={id} className="col-12 text-center">
                  <Button
                    variant="link"
                    className={
                      id == displayTab && Object.keys(expVarList).length
                        ? 'secondary-navlinks px-3 py-1 d-inline-block border-0 bg-association text-white rounded-0'
                        : 'secondary-navlinks px-3 py-1 d-inline-block border-0 rounded-0'
                    }
                    style={{
                      textDecoration: 'none',
                      fontSize: '12pt',
                      color: '#b83d47',
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
            {/* <Row>
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
            </Row> */}
            <Row>
              <Col lg="auto" className="w-100">
                {source == 'public' ? <PublicForm /> : <></>}
                {/* <UserForm /> */}
              </Col>
            </Row>
          </div>
          <hr className="d-lg-none" style={{ opacity: 0 }}></hr>
        </SidebarPanel>
        <MainPanel>
          <div style={{ minHeight: '500px' }}>
            {tabs.filter((tab) => tab.id == displayTab)[0].component}
          </div>
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
