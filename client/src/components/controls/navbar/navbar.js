import React from 'react';
import { NavLink } from 'react-router-dom';
import { Navbar, Container, Nav } from 'react-bootstrap';

import './navbar.scss';

export function NavbarCustom({ links }) {
  return (
    <div>
      <Navbar
        collapseOnSelect
        className="shadow-sm bg-nav p-0 border-0 d-none d-lg-block"
        expand="lg"
      >
        <div className="container">
          <Navbar.Toggle aria-controls="basic-navbar-nav" />
          <Navbar.Collapse id="basic-navbar-nav">
            <div className="mx-auto">
              {[{ route: '/', title: 'Home', exact: true }]
                .concat(links)
                .filter((link) => link.title)
                .sort((a, b) => a.navIndex - b.navIndex)
                .map(({ route, action, title, exact }, index) => (
                  <div
                    data-testid="Navbar"
                    className="d-inline-block"
                    key={title}
                  >
                    <NavLink
                      data-testid={`Navbar-NavLink-${index}`}
                      id={title + '-navbar'}
                      // key={title}
                      className="navlinks py-2 px-3 d-inline-block"
                      activeClassName="active-navlinks"
                      exact={exact}
                      to={route}
                    >
                      {title}
                    </NavLink>
                  </div>
                ))}
            </div>
          </Navbar.Collapse>
        </div>
      </Navbar>

      {/* Mobile view*/}
      <Navbar
        collapseOnSelect
        className="shadow-sm bg-nav border-0 d-block d-lg-none ml-3 mr-auto"
        expand="lg"
      >
        <div className="mx-auto">
          <Navbar.Toggle aria-controls="basic-navbar-nav" />
          <Navbar.Collapse id="basic-navbar-nav">
            <div className="mx-auto">
              {[{ route: '/', title: 'Home', exact: true }]
                .concat(links)
                .filter((link) => link.title)
                .sort((a, b) => a.navIndex - b.navIndex)
                .map(({ route, action, title, exact }, index) => (
                  <div data-testid="Navbar" className="" key={title}>
                    <NavLink
                      data-testid={`Navbar-NavLink-${index}`}
                      id={title + '-navbar'}
                      // key={title}
                      className="navlinks py-2 px-4 d-inline-block"
                      activeClassName="active-navlinks"
                      exact={exact}
                      to={route}
                    >
                      {title}
                    </NavLink>
                  </div>
                ))}
            </div>
          </Navbar.Collapse>
        </div>
      </Navbar>
    </div>
  );
}
