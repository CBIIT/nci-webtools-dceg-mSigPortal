import React, { useState } from 'react';
import { Accordion, Card } from 'react-bootstrap';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faMinus, faPlus } from '@fortawesome/free-solid-svg-icons';
import './accordions.scss';

const { Header, Body } = Card;
const { Toggle, Collapse } = Accordion;

export default function ({ components }) {
  const [toggleState, setToggle] = useState({});

  return (
    <div className="accordions">
      {components.map(({ title, component }, index) => {
        return (
          <Accordion defaultActiveKey={index + ''} key={index + ''}>
            <Card>
              <Toggle
                className="font-weight-bold"
                as={Header}
                eventKey={index + ''}
                onClick={() =>
                  setToggle({
                    ...toggleState,
                    [index + '']: !toggleState[index + ''],
                  })
                }
              >
                {toggleState[index + ''] ? (
                  <FontAwesomeIcon icon={faPlus} />
                ) : (
                  <FontAwesomeIcon icon={faMinus} />
                )}{' '}
                {title}
              </Toggle>
              <Collapse eventKey={index + ''}>
                <Body>{component}</Body>
              </Collapse>
            </Card>
          </Accordion>
        );
      })}
    </div>
  );
}
