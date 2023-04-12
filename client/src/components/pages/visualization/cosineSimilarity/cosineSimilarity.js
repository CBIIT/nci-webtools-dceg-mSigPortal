import { useState } from 'react';
import { Nav, Tab } from 'react-bootstrap';
import CsWithin from './cosineSimilarity-within';
import CsReference from './cosineSimilarity-reference';
import CsPublic from './cosineSimilarity-public';
import { NavHashLink } from 'react-router-hash-link';

export default function CosineSimilarity({ state }) {
  const { source } = state;
  const [display, setDisplay] = useState('within');

  return (
    <div>
      <p className="bg-white border rounded p-3 mb-1">
        Below you can explore cosine similarity between sample profiles (CS
        Between Samples), cosine similarity between sample profiles and
        reference signatures (CS to Reference Signatures), or, if using your own
        data, cosine similarity between profiles from your input data and
        profiles from public data (CS to Public Data). Simply use the dropdown
        menus to select a [Profile Type], [Matrix Size], or [Reference Signature
        Set]. Click here to learn more about cosine similarity. Click{' '}
        <NavHashLink to="/faq#cosine-similarity">here</NavHashLink> to learn
        more about cosine similarity.
      </p>
      <Tab.Container
        transition={false}
        className="mt-2"
        defaultActiveKey={display}
        activeKey={display}
        onSelect={(tab) => setDisplay(tab)}
      >
        <Nav variant="tabs">
          <Nav.Item>
            <Nav.Link eventKey="within" as="button" className="outline-none">
              <strong>CS Between Samples</strong>
            </Nav.Link>
          </Nav.Item>
          <Nav.Item>
            <Nav.Link eventKey="reference" as="button" className="outline-none">
              <strong>CS to Reference Signatures</strong>
            </Nav.Link>
          </Nav.Item>
          {source == 'user' && (
            <Nav.Item>
              <Nav.Link eventKey="public" as="button" className="outline-none">
                <strong>CS to Public Data</strong>
              </Nav.Link>
            </Nav.Item>
          )}
        </Nav>
        <Tab.Content
          className={`bg-white tab-pane-bordered rounded-0 d-block`}
          style={{ overflowX: 'auto' }}
        >
          <Tab.Pane key="within" eventKey="within" className="border-0">
            <CsWithin state={state} />
          </Tab.Pane>
          <Tab.Pane key="reference" eventKey="reference" className="border-0">
            <CsReference state={state} />
          </Tab.Pane>
          {source == 'user' && (
            <Tab.Pane key="public" eventKey="public" className="border-0">
              <CsPublic state={state} />
            </Tab.Pane>
          )}
        </Tab.Content>
      </Tab.Container>
    </div>
  );
}
