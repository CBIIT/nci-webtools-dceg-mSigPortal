import { Button, Nav } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import InputForm from './inputForm/inputForm';
import Instructions from './instructions';
// import Download from './download';

import { actions as extractionActions } from '../../../services/store/extraction';
import { actions as modalActions } from '../../../services/store/modal';
import {
  SidebarContainer,
  SidebarPanel,
  MainPanel,
} from '../../controls/sidebar-container/sidebar-container';

const actions = { ...extractionActions, ...modalActions };

export default function Extraction({ match }) {
  const dispatch = useDispatch();
  const mergeState = (state) =>
    dispatch(actions.mergeExtraction({ main: state }));

  const { displayTab, openSidebar, submitted } = useSelector(
    (state) => state.extraction.main
  );

  const tabs = [
    {
      component: <Instructions loading={false} />,
      id: 'instructions',
      name: 'Instructions',
    },
  ];

  return (
    <div className="position-relative">
      <div className="mx-3">
        <div className="mx-3 bg-white border border-top-0">
          {/* for desktops and tablets */}
          <div className="d-none d-md-block">
            <Nav defaultActiveKey="instructions">
              {tabs.map(({ name, id }) => {
                if (name)
                  return (
                    <div key={id} className="d-inline-block">
                      <Button
                        variant="link"
                        className={`secondary-navlinks px-3 py-1 d-inline-block border-0 ${
                          id == displayTab ? 'active-secondary-navlinks' : ''
                        }`}
                        active={id == displayTab && submitted}
                        disabled={id != 'instructions' && !submitted}
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
                    </div>
                  );
              })}
            </Nav>
          </div>
          {/* for mobile devices */}
          <div className="row d-md-none">
            <Nav defaultActiveKey="instructions">
              {tabs.map(({ name, id }) => {
                if (name)
                  return (
                    <div key={id} className="col-12 text-center">
                      <Button
                        variant="link"
                        className={
                          'secondary-navlinks px-3 py-1 d-inline-block border-0'
                        }
                        active={id == displayTab && submitted}
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
            <InputForm />
          </div>
        </SidebarPanel>
        <MainPanel>
          {tabs.filter((tab) => tab.id == displayTab)[0].component}
        </MainPanel>
      </SidebarContainer>
    </div>
  );
}
