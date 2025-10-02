import React, { useState } from 'react';
import { Button, Nav } from 'react-bootstrap';
import Instructions from './instructions';
import Status from './status';
import TargetedSequencing from './targetedSequencing';

export default function Refitting() {
  const [displayTab, setDisplayTab] = useState('instructions');

  const tabs = [
    {
      id: 'instructions',
      name: 'Instructions',
      disabled: false,
    },
    {
      id: 'status',
      name: 'Status',
      disabled: false,
    },
    {
      id: 'targetedSequencing',
      name: 'Targeted Sequencing',
      disabled: false, // You might want to make this conditional based on analysis status
    },
  ];

  return (
    <div className="position-relative">
      <div className="mx-3">
        <div className="mx-3 bg-white border border-top-0">
          {/* for desktops and tablets */}
          <div className="d-none d-md-block">
            <Nav defaultActiveKey="instructions">
              {tabs.map(({ name, id, disabled }, i) => (
                <div key={name} className="d-inline-block">
                  <Button
                    variant="link"
                    className={`secondary-navlinks px-3 py-1 d-inline-block border-0 text-refitting rounded-0 ${
                      id === displayTab
                        ? 'bg-refitting text-white'
                        : 'text-refitting'
                    }`}
                    disabled={disabled}
                    style={{
                      textDecoration: 'none',
                      fontSize: '12pt',
                      color: '#689f39',
                      fontWeight: '500',
                    }}
                    onClick={() => setDisplayTab(id)}
                  >
                    {name}
                  </Button>
                </div>
              ))}
            </Nav>
          </div>

          {/* for mobile devices */}
          <div className="row d-md-none">
            <Nav defaultActiveKey="instructions">
              {tabs.map(({ name, id, disabled }, i) => {
                if (name)
                  return (
                    <div key={name} className="col-12 text-center">
                      <Button
                        variant="link"
                        className={
                          id === displayTab
                            ? 'secondary-navlinks px-3 py-1 d-inline-block border-0 bg-refitting text-white rounded-0'
                            : 'secondary-navlinks px-3 py-1 d-inline-block border-0 rounded-0'
                        }
                        disabled={disabled}
                        style={{
                          textDecoration: 'none',
                          fontSize: '12pt',
                          color: '#689f39',
                          fontWeight: '500',
                        }}
                        onClick={() => setDisplayTab(id)}
                      >
                        {name}
                      </Button>
                      <div className="d-md-none w-100"></div>
                    </div>
                  );
              })}
            </Nav>
          </div>
        </div>
      </div>

      <div className="m-3">
        <div className={displayTab === 'instructions' ? 'd-block' : 'd-none'}>
          <Instructions />
        </div>
        <div className={displayTab === 'status' ? 'd-block' : 'd-none'}>
          <Status />
        </div>
        <div className={displayTab === 'targetedSequencing' ? 'd-block' : 'd-none'}>
          <TargetedSequencing />
        </div>
      </div>
    </div>
  );
}
