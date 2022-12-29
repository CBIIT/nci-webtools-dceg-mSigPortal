import React from 'react';
import { NavLink } from 'react-router-dom';
import { Navbar, Container, Nav } from 'react-bootstrap';

import './navbar.scss';

export function NavbarCustom({ links }) {
  return (
    <div>
      <div className="shadow-sm bg-nav border-0">
        <div className="d-none d-md-flex justify-content-center">
          {[{ route: '/', title: 'Home', exact: true }]
            .concat(links)
            .filter((link) => link.title)
            .sort((a, b) => a.navIndex - b.navIndex)
            .map(({ route, action, title, exact }, index) => (
              <div data-testid="Navbar" className="d-inline-block" key={title}>
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
          {/* <pre>{JSON.stringify(links)}</pre> */}
        </div>

        {/*Mobile View*/}
        <div className="container d-flex flex-column d-md-none" align="center">
          {[{ route: '/', title: 'Home', exact: true }]
            .concat(links)
            .filter((link) => link.title)
            .sort((a, b) => a.navIndex - b.navIndex)
            .map(({ route, action, title, exact, dropdown }, index) => {
              if (dropdown) {
                return dropdown.map(({ name, path }) => (
                  <div data-testid="Navbar" className="d-block" key={name}>
                    <NavLink
                      data-testid={`Navbar-NavLink-${index}`}
                      id={'nav' + name}
                      // key={title}
                      className="navlinks py-2 px-4 d-inline-block"
                      activeClassName="active-navlinks"
                      exact={exact}
                      to={`/catalog/${path}`}
                    >
                      {name}
                    </NavLink>
                  </div>
                ));
              } else {
                return (
                  <div data-testid="Navbar" className="d-block" key={title}>
                    <NavLink
                      data-testid={`Navbar-NavLink-${index}`}
                      id={'nav' + title}
                      // key={title}
                      className="navlinks py-2 px-4 d-inline-block"
                      activeClassName="active-navlinks"
                      exact={exact}
                      to={route}
                    >
                      {title}
                    </NavLink>
                  </div>
                );
              }
            })}
        </div>
      </div>
    </div>
  );
}