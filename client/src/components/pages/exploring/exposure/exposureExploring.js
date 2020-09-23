import React from 'react';
import { Row, Col, Accordion, Card } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import Tumor from './tumor';
import Activity from './activity';
import Decomposition from './decomposition';
import Landscape from './landscape';
import Prevalence from './prevalence';
import { dispatchError, dispatchExploring } from '../../../../services/store';

const { Header, Body } = Card;
const { Toggle, Collapse } = Accordion;

export default function ExposureExploring() {
  const rootURL = window.location.pathname;
  const { displayTab, projectID, exposureAccordion } = useSelector(
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

  const sections = [
    {
      component: <Tumor submitR={(fn, args) => submitR(fn, args)} />,
      id: 'tumor',
      title: 'Tumor Mutational Burden',
    },
    {
      component: <Activity submitR={(fn, args) => submitR(fn, args)} />,
      id: 'activity',
      title: 'Mutational Signature Activity',
    },
    {
      component: (
        <Decomposition
          downloadResults={(path) => downloadResults(path)}
          submitR={(fn, args) => submitR(fn, args)}
        />
      ),
      id: 'decomposition',
      title: 'Evaluating the Performance of Mutational Signature Decomposition',
    },
    {
      component: <Landscape submitR={(fn, args) => submitR(fn, args)} />,
      id: 'landscape',
      title: 'Landscape of Mutational Signature Activity',
    },
    {
      component: <Prevalence submitR={(fn, args) => submitR(fn, args)} />,
      id: 'prevalence',
      title: 'Prevalence of Mutational Signature',
    },
  ];
  return (
    <div className="position-relative">
      {sections.map(({ component, id, title }) => {
        return (
          <Accordion activeKey={exposureAccordion[id]} key={id}>
            <Card>
              <Toggle
                className="font-weight-bold"
                as={Header}
                eventKey={exposureAccordion[id]}
                onClick={() =>
                  dispatchExploring({
                    exposureAccordion: {
                      ...exposureAccordion,
                      [id]: !exposureAccordion[id],
                    },
                  })
                }
              >
                {exposureAccordion[id] == true ? (
                  <FontAwesomeIcon icon={faMinus} />
                ) : (
                  <FontAwesomeIcon icon={faPlus} />
                )}{' '}
                {title}
              </Toggle>
              <Collapse eventKey={exposureAccordion[id]}>
                <Body>{component}</Body>
              </Collapse>
            </Card>
          </Accordion>
        );
      })}
    </div>
  );
}
