import React from 'react';
import { Row, Col, Accordion, Card } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import ReferenceSignatures from './referenceSignatures';
import MutationalSignatureProfile from './mutationalSignatureProfile';
import CosineSimilarity from './cosineSimilarity';
import MutationalSignatureComparison from './mutationalSignatureComparison';
import { dispatchError, dispatchExploring } from '../../../../services/store';

const { Header, Body } = Card;
const { Toggle, Collapse } = Accordion;

export default function SignatureExploring() {
  const rootURL = window.location.pathname;

  const {
    displayTab,
    projectID,
    signatureAccordion,
    exposureAccordion,
  } = useSelector((state) => state.exploring);

  function submitR(fn, args) {
    return fetch(`${rootURL}exploringR`, {
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
    <div className="position-relative">
      {sections.map(({ component, id, title }) => {
        return (
          <Accordion activeKey={signatureAccordion[id]} key={id}>
            <Card>
              <Toggle
                className="font-weight-bold"
                as={Header}
                eventKey={signatureAccordion[id]}
                onClick={() =>
                  dispatchExploring({
                    signatureAccordion: {
                      ...signatureAccordion,
                      [id]: !signatureAccordion[id],
                    },
                  })
                }
              >
                {signatureAccordion[id] == true ? (
                  <FontAwesomeIcon icon={faMinus} />
                ) : (
                  <FontAwesomeIcon icon={faPlus} />
                )}{' '}
                {title}
              </Toggle>
              <Collapse eventKey={signatureAccordion[id]}>
                <Body>{component}</Body>
              </Collapse>
            </Card>
          </Accordion>
        );
      })}
    </div>
  );
}
