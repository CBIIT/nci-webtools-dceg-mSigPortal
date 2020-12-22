import React, { useEffect } from 'react';
import { Card, Nav } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import ReferenceSignatures from './referenceSignatures';
import MutationalSignatureProfile from './mutationalSignatureProfile';
import CosineSimilarity from './cosineSimilarity';
import MutationalSignatureComparison from './mutationalSignatureComparison';
import { dispatchExploring } from '../../../../services/store';
import Accordions from '../../../controls/accordions/accordions';

const { Header, Body } = Card;

export default function SignatureExploring() {
  const { projectID } = useSelector((state) => state.exploring);

  useEffect(() => {
    dispatchExploring({ displayTab: 'signature' });
  }, []);

  function submitR(fn, args) {
    return fetch(`api/exploringR`, {
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

  const sections = [
    {
      component: (
        <ReferenceSignatures submitR={(fn, args) => submitR(fn, args)} />
      ),
      id: 'referenceSignatures',
      title: 'Current Reference Signatures in mSigPortal',
    },
    {
      component: (
        <MutationalSignatureProfile submitR={(fn, args) => submitR(fn, args)} />
      ),
      id: 'mutationalSignatureProfile',
      title: 'Mutational Signature Profile',
    },
    {
      component: <CosineSimilarity submitR={(fn, args) => submitR(fn, args)} />,
      id: 'cosineSimilarity',
      title: 'Cosine Similarity Among Mutational Signatures',
    },
    {
      component: (
        <MutationalSignatureComparison
          submitR={(fn, args) => submitR(fn, args)}
        />
      ),
      id: 'mutationalSignatureComparison',
      title: 'Mutational Signatures Comparisons',
    },
  ];

  return (
    <Card>
      <Header>
        <Nav variant="pills" defaultActiveKey="#exploring/signature">
          <Nav.Item>
            <Nav.Link href="#exploring/signature">Signature Exploring</Nav.Link>
          </Nav.Item>
        </Nav>
      </Header>
      <Body>
        <Accordions components={sections} />
      </Body>
    </Card>
  );
}
