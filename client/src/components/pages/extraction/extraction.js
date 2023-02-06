import { Button, Nav } from 'react-bootstrap';
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

const actions = { ...extractionActions, ...modalActions };

export default function Extraction() {
  const dispatch = useDispatch();
  const mergeState = (state) => dispatch(actions.mergeExtraction(state));

  const id = useParams().id || false;
  const { tabIndex, openSidebar, submitted } = useSelector(
    (state) => state.extraction
  );

  const tabs = [
    {
      name: 'Instructions',
      component: <Instructions loading={false} />,
    },
    {
      name: 'TMB',
      component: <TMB />,
    },
  ];

  return (
    <div className="position-relative">
      <div className="mx-3">
        <div className="mx-3 bg-white border border-top-0">
          {/* for desktops and tablets */}
          <div className="d-none d-md-block">
            <Nav defaultActiveKey="instructions">
              {tabs.map(({ name }, i) => (
                <div key={name} className="d-inline-block">
                  <Button
                    variant="link"
                    className={`secondary-navlinks px-3 py-1 d-inline-block border-0 rounded-0 ${
                      i == tabIndex ? 'bg-extraction text-white' : ''
                    }`}
                    disabled={name != 'Instructions' && !submitted}
                    style={{
                      textDecoration: 'none',
                      fontSize: '12pt',
                      color: '#42688b',
                      fontWeight: '500',
                    }}
                    onClick={() => mergeState({ displayTab: i })}
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
              {tabs.map(({ name }, i) => {
                if (name)
                  return (
                    <div key={name} className="col-12 text-center">
                      <Button
                        variant="link"
                        className={
                          'secondary-navlinks px-3 py-1 d-inline-block border-0 bg-extraction text-white'
                        }
                        style={{
                          textDecoration: 'none',
                          fontSize: '12pt',
                          color: '#42688b',
                          fontWeight: '500',
                        }}
                        onClick={() => mergeState({ displayTab: i })}
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
        <MainPanel>{tabs[tabIndex].component}</MainPanel>
      </SidebarContainer>
    </div>
  );
}
