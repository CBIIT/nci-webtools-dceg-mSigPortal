import React, { useEffect } from 'react';
import { Tab, Nav } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import ReferenceSignatures from './referenceSignatures';
import MutationalSignatureProfile from './mutationalSignatureProfile';
import CosineSimilarity from './cosineSimilarity';
import MutationalSignatureComparison from './mutationalSignatureComparison';
import {
  actions,
  getInitialState,
} from '../../../../services/store/exploration';

const { Container, Content, Pane } = Tab;
const { Item, Link } = Nav;

export default function Signature() {
  const dispatch = useDispatch();
  const exploration = useSelector((state) => state.exploration);
  const mergeExploration = (state) =>
    dispatch(actions.mergeExploration({ exploration: state }));

  const { projectID, signatureDisplay } = exploration.exploration;

  useEffect(() => {
    mergeExploration({ displayTab: 'signature' });
  }, []);

  function submitR(fn, args) {
    return fetch(`api/explorationCalc`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        fn: fn,
        args: args,
        projectID: projectID,
      }),
    });
  }

  const tabs = [
    {
      component: (
        <ReferenceSignatures submitR={(fn, args) => submitR(fn, args)} />
      ),
      key: 'referenceSignatures',
      title: 'Current Reference Signatures in mSigPortal',
    },
    {
      component: (
        <MutationalSignatureProfile submitR={(fn, args) => submitR(fn, args)} />
      ),
      key: 'mutationalSignatureProfile',
      title: 'Mutational Signature Profile',
    },
    {
      component: <CosineSimilarity submitR={(fn, args) => submitR(fn, args)} />,
      key: 'cosineSimilarity',
      title: 'Cosine Similarity Among Mutational Signatures',
    },
    {
      component: (
        <MutationalSignatureComparison
          submitR={(fn, args) => submitR(fn, args)}
        />
      ),
      key: 'mutationalSignatureComparison',
      title: 'Mutational Signatures Comparisons',
    },
  ];

  return (
    <div>
      <Container
        transition={false}
        className="mt-2"
        defaultActiveKey={signatureDisplay}
        activeKey={signatureDisplay}
        onSelect={(tab) => mergeExploration({ signatureDisplay: tab })}
      >
        <Nav variant="tabs">
          {tabs.map(({ key, title }) => (
            <Item key={key}>
              <Link eventKey={key} as="button" className="outline-none">
                <strong>{title}</strong>
              </Link>
            </Item>
          ))}
        </Nav>
        <Content
          className={`bg-white tab-pane-bordered rounded-0 d-block`}
          style={{ overflowX: 'auto' }}
        >
          {tabs.map(({ key, component }) => (
            <Pane key={key} eventKey={key} className="border-0">
              {component}
            </Pane>
          ))}
        </Content>
      </Container>
    </div>
  );
}
