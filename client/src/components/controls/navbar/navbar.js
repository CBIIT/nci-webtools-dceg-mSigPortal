import React, { useState } from 'react';
import { NavLink, useLocation } from 'react-router-dom';
import { NavDropdown } from 'react-bootstrap';
import './navbar.scss';

export function Navbar({ links }) {
  const [dropdownStatus, setDropdown] = useState(false);
  const location = useLocation();

  return (
    <div className="bg-dark text-white gradient shadow-sm">
      <div className="container d-none d-md-flex justify-content-center">
        {[{ route: '/', title: 'Home', exact: true }]
          .concat(links)
          .filter((link) => link.title)
          .sort((a, b) => a.navIndex - b.navIndex)
          .map(({ route, action, title, exact, dropdown }, index) => {
            if (dropdown) {
              return (
                <NavDropdown
                  key={title}
                  id={title}
                  title={title}
                  active={location.pathname.indexOf(route) > -1}
                  onMouseEnter={() => setDropdown(true)}
                  onMouseLeave={() => setDropdown(false)}
                  show={dropdownStatus}
                >
                  {dropdown.map(({ name, path }) => (
                    <NavLink
                      activeClassName="active-navlinks"
                      key={`${title}/${path}`}
                      className="dropdown-item"
                      exact={true}
                      to={`/exploring/${path}`}
                    >
                      {name}
                    </NavLink>
                  ))}
                </NavDropdown>
              );
            } else {
              return (
                <div
                  data-testid="Navbar"
                  className="d-inline-block"
                  key={title}
                >
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
              );
            }
          })}
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
