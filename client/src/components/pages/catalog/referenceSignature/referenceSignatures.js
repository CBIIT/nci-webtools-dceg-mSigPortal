import { Tab, Nav } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import RsInMsigportal from './rsInMsigportal/rsInMsigportal-plot';
import RsProfile from './rsProfile/rsProfile';
import CosineSimilarity from './cosineSimilarity/cosineSimilarity';
import RsComparison from './rsComparison/rsComparison';
import Download from './download';
import { actions } from '../../../../services/store/catalog';

export default function ReferenceSignature() {
  const dispatch = useDispatch();
  const { displayTab } = useSelector(
    (state) => state.catalog.referenceSignature
  );
  const mergeState = (state) =>
    dispatch(actions.mergeCatalog({ referenceSignature: state }));

  const tabs = [
    {
      component: <RsInMsigportal />,
      key: 'overview',
      title: 'RS In mSigPortal',
    },
    {
      component: <RsProfile />,
      key: 'RSProfile',
      title: 'RS Profile',
    },
    {
      component: <CosineSimilarity />,
      key: 'cosineSimilarity',
      title: 'RS Cosine Similarity',
    },
    {
      component: <RsComparison />,
      key: 'SignatureComparison',
      title: 'RS Comparison',
    },
    {
      component: <Download />,
      key: 'download',
      title: 'Download',
    },
  ];

  return (
    <div style={{ minHeight: '500px' }}>
      <Tab.Container
        transition={false}
        className="mt-2"
        activeKey={displayTab}
        onSelect={(tab) => mergeState({ displayTab: tab })}
      >
        <Nav variant="tabs">
          {tabs.map(({ key, title }) => (
            <Nav.Item key={key}>
              <Nav.Link
                eventKey={key}
                as="button"
                className="outline-none catalog"
              >
                <strong>{title}</strong>
              </Nav.Link>
            </Nav.Item>
          ))}
        </Nav>
        <Tab.Content className="bg-white tab-pane-bordered rounded-0 d-block">
          {tabs.map(({ key, component }) => (
            <Tab.Pane key={key} eventKey={key} className="border-0">
              {component}
            </Tab.Pane>
          ))}
        </Tab.Content>
      </Tab.Container>
    </div>
  );
}
