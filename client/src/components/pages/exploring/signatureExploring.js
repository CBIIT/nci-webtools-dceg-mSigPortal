import React from 'react';
import { Row, Col, Accordion, Card } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import ReferenceSignatures from './referenceSignatures';
import MutationalSignatureProfile from './mutationalSignatureProfile';
import CosineSimilarity from './cosineSimilarity';
import { dispatchError, dispatchExploring } from '../../../services/store';

const { Header, Body } = Card;
const { Toggle, Collapse } = Accordion;

export default function SignatureExploring({ getReferenceSignatureData }) {
  const rootURL = window.location.pathname;

  const { displayTab, signatureTab, exposureTab, projectID } = useSelector(
    (state) => state.exploring
  );

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

  //   download text results files
  async function downloadResults(txtPath) {
    try {
      const response = await fetch(`${rootURL}downloadPlotData`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ path: txtPath }),
      });
      if (!response.ok) {
        const { msg } = await response.json();
        dispatchError(msg);
      } else {
        const file = await response.blob();
        const objectURL = URL.createObjectURL(file);
        const tempLink = document.createElement('a');

        tempLink.href = objectURL;
        tempLink.download = txtPath.split('/').slice(-1)[0];
        document.body.appendChild(tempLink);
        tempLink.click();
        document.body.removeChild(tempLink);
        URL.revokeObjectURL(objectURL);
      }
    } catch (err) {
      dispatchError(err);
    }
  }

  function getRefSigOptions(profileType) {
    return fetch(`${rootURL}visualizeR/getReferenceSignatureSets`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ profileType: profileType }),
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
        <MutationalSignatureProfile
          submitR={(fn, args) => submitR(fn, args)}
          getReferenceSignatureData={(c, f) => getReferenceSignatureData(c, f)}
        />
      ),
      id: 'mutationalSignatureProfile',
      title: 'Mutational Signature Profile',
    },
    {
      component: (
        <CosineSimilarity
          getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
          downloadResults={(path) => downloadResults(path)}
          submitR={(fn, args) => submitR(fn, args)}
          getReferenceSignatureData={(c, f) => getReferenceSignatureData(c, f)}
        />
      ),
      id: 'cosineSimilarity',
      title: 'Cosine Similarity Among Mutational Signatures',
    },
  ];
  return (
    <div className="position-relative">
      {sections.map(({ component, id, title }) => {
        return (
          <Accordion defaultActiveKey={id} key={id}>
            <Card>
              <Toggle
                className="font-weight-bold"
                as={Header}
                eventKey={id}
                onClick={() =>
                  dispatchExploring({
                    signatureTab: {
                      ...signatureTab,
                      [id]: !signatureTab[id],
                    },
                  })
                }
              >
                {signatureTab[id] == true ? (
                  <FontAwesomeIcon icon={faMinus} />
                ) : (
                  <FontAwesomeIcon icon={faPlus} />
                )}{' '}
                {title}
              </Toggle>
              <Collapse eventKey={id}>
                <Body>{component}</Body>
              </Collapse>
            </Card>
          </Accordion>
        );
      })}
    </div>
  );
}
