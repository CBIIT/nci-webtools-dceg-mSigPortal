import { Nav, Tab } from 'react-bootstrap';
import PcWithin from './profileComparison-within';
import PcReference from './profileComparison-reference';
import PcPublic from './profileComparison-public';
import { useSelector, useDispatch } from 'react-redux';
import { actions } from '../../../../services/store/visualization';

export default function MutationalProfiles(props) {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.visualization);
  const mergeState = (state) =>
    dispatch(actions.mergeVisualization({ profileComparison: state }));

  const { source } = store.main;
  const { display } = store.profileComparison;

  return (
    <div>
      <p className="bg-white border rounded p-3 mb-1">
        Below you can perform mutational profile comparison analyses between
        samples (PC Between Samples), between a sample and a signature from a
        selected reference signature set (PC to Reference Signatures), or
        between user data and public data (PC to Public Data).
      </p>
      <Tab.Container
        transition={false}
        className="mt-2"
        defaultActiveKey={display}
        activeKey={display}
        onSelect={(tab) => mergeState({ display: tab })}
      >
        <Nav variant="tabs">
          <Nav.Item>
            <Nav.Link eventKey="within" as="button" className="outline-none">
              <strong>PC Between Samples</strong>
            </Nav.Link>
          </Nav.Item>
          <Nav.Item>
            <Nav.Link eventKey="reference" as="button" className="outline-none">
              <strong>PC to Reference Signatures</strong>
            </Nav.Link>
          </Nav.Item>
          {source == 'user' && (
            <Nav.Item>
              <Nav.Link eventKey="public" as="button" className="outline-none">
                <strong>PC to Public Data</strong>
              </Nav.Link>
            </Nav.Item>
          )}
        </Nav>
        <Tab.Content
          className={`bg-white tab-pane-bordered rounded-0 d-block`}
          style={{ overflowX: 'auto' }}
        >
          <Tab.Pane key="within" eventKey="within" className="border-0">
            <PcWithin />
          </Tab.Pane>
          <Tab.Pane key="reference" eventKey="reference" className="border-0">
            <PcReference />
          </Tab.Pane>
          {source == 'user' && (
            <Tab.Pane key="public" eventKey="public" className="border-0">
              <PcPublic />
            </Tab.Pane>
          )}
        </Tab.Content>
      </Tab.Container>
    </div>
  );
}
