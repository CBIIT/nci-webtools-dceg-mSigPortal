import { Nav, Tab } from 'react-bootstrap';
import PcaWithin from './pca-within';
import PcaPublic from './pca-public';
import { useSelector, useDispatch } from 'react-redux';
import Description from '../../../controls/description/description';
import { actions } from '../../../../services/store/visualization';

export default function PCA() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.visualization);
  const mergeState = (state) =>
    dispatch(actions.mergeVisualization({ pca: state }));

  const { source } = store.main;
  const { display } = store.pca;

  return (
    <div>
      <div className="bg-white border rounded p-3 mb-1">
        <Description
          less="Below you can conduct PCA analysis between samples, or a PCA with Public Data (for user input data only)."
          more="PCA stands for Principal Component Analysis, which helps to explain the variation found in the data through the establishment of different principal components. Each principal component can also be used to compare with known mutational signatures."
        />
      </div>
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
              <strong>PCA Between Samples</strong>
            </Nav.Link>
          </Nav.Item>
          {source == 'user' && (
            <Nav.Item>
              <Nav.Link eventKey="public" as="button" className="outline-none">
                <strong>PCA to Public Data</strong>
              </Nav.Link>
            </Nav.Item>
          )}
        </Nav>
        <Tab.Content
          className={`bg-white tab-pane-bordered rounded-0 d-block`}
          style={{ overflowX: 'auto' }}
        >
          <Tab.Pane key="within" eventKey="within" className="border-0">
            <PcaWithin />
          </Tab.Pane>
          {source == 'user' && (
            <Tab.Pane key="public" eventKey="public" className="border-0">
              <PcaPublic />
            </Tab.Pane>
          )}
        </Tab.Content>
      </Tab.Container>
    </div>
  );
}
