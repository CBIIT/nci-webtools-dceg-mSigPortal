import React from 'react';
import { Nav } from 'react-bootstrap';
import { Route, Redirect, NavLink } from 'react-router-dom';
import { useSelector, useDispatch } from 'react-redux';
import ReferenceSignature from './referenceSignature/referenceSignatures';
import Etiology from './etiology/etiology';
import { actions } from '../../../services/store/catalog';

export default function Explore() {
  const dispatch = useDispatch();
  const mergeState = (state) => dispatch(actions.mergeCatalog({ main: state }));

  const store = useSelector((state) => state.catalog);
  const { displayTab } = store.main;

  const tabs = [
    { name: 'Signature Catalog', id: 'etiology' },
    { name: 'Reference Signature', id: 'referenceSignature' },
  ];

  function handleTabChange(tab) {
    mergeState({ displayTab: tab });
  }

  return (
    <div className="px-0">
      <div className="">
        <div className="d-flex justify-content-center border border-catalog border-top-0 border-left-0 border-right-0">
          {/* for desktops and tablets */}
          <div className="d-none d-md-block">
            <Nav activeKey={displayTab}>
              {tabs.map(({ name, id }) => (
                <Nav.Item key={id} className="d-inline-block">
                  <NavLink
                    className="secondary-navlinks px-3 py-1 d-inline-block catalogNav text-catalog"
                    activeClassName="bg-catalog text-white"
                    style={{
                      textDecoration: 'none',
                      fontSize: '12pt',
                      //color: 'black',
                      fontWeight: '500',
                    }}
                    exact={true}
                    to={`/catalog/${id}`}
                    onClick={() => handleTabChange(id)}
                  >
                    {name}
                  </NavLink>
                  <div className="d-md-none w-100"></div>
                </Nav.Item>
              ))}
            </Nav>
          </div>
          {/* for mobile devices */}
          <div className="row d-md-none">
            <Nav activeKey={displayTab}>
              {tabs.map(({ name, id }) => (
                <Nav.Item key={id} className="col-12 text-center">
                  <NavLink
                    className="secondary-navlinks px-3 py-1 d-inline-block text-catalog"
                    activeClassName="bg-catalog text-white"
                    style={{
                      textDecoration: 'none',
                      fontSize: '12pt',
                      //color: 'black',
                      fontWeight: '500',
                    }}
                    exact={true}
                    to={`/catalog/${id}`}
                    onClick={() => handleTabChange(id)}
                  >
                    {name}
                  </NavLink>
                  <div className="d-md-none w-100"></div>
                </Nav.Item>
              ))}
            </Nav>
          </div>
        </div>
        <div className="mx-3 my-3">
          <Route
            exact
            path={`/catalog`}
            render={() => <Redirect to={`/catalog/${displayTab}`} />}
          />
          <Route path="/catalog/etiology" component={Etiology} />
          <Route
            path="/catalog/referenceSignature"
            component={ReferenceSignature}
          />
        </div>
      </div>
    </div>
  );
}
