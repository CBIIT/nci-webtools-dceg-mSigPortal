import React from 'react';
import { NavLink } from 'react-router-dom';
import './navbar.scss';

export function Navbar({ links }) {
  return (
    <div className="bg-dark text-white gradient shadow-sm">
      <div className="container d-none d-md-flex justify-content-center">
        {[{ route: '/', title: 'Home', exact: true }]
          .concat(links)
          .filter((link) => link.title)
          .sort((a, b) => a.navIndex - b.navIndex)
          .map(({ route, action, title, exact }, index) => (
            <div data-testid="Navbar" className="d-inline-block" key={title}>
              <NavLink
                data-testid={`Navbar-NavLink-${index}`}
                id={title}
                // key={title}
                className="navlinks py-2 px-4 d-inline-block"
                activeClassName="active-navlinks"
                exact={exact}
                to={route}
              >
                {title}
              </NavLink>
              <div className="d-lg-none w-100"></div>
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
          .map(({ route, action, title, exact }, index) => (
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
          ))}
      </div>
    </div>
  );
}
