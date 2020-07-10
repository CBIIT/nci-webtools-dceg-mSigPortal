import React, { useState, useEffect } from 'react';
import { Link } from 'react-router-dom';
import Card from 'react-bootstrap/Card';
import { CardDeck, Button } from 'react-bootstrap';
import './home.scss';

export default function Home({ links }) {
  const colors = ['#fc8701', '#2c71dd', '#689f39', '#84368d'];

  function cardRow(links, colors) {
    return (
      <CardDeck>
        {links.map(
          (
            { exact, route, action, title, cardTitle, cardText, image },
            index
          ) => (
            <>
              <Card
                key={title}
                id={title}
                className="mb-5 align-self-center"
                style={{
                  width: '18rem',
                  // height: '18rem',
                  justifyContent: 'center',
                  alignItems: 'center',
                  border: '1px solid #DADBE6',
                  backgroundColor: 'white',
                  // minHeight: '280px'
                  // borderRadius: '10px'
                }}
              >
                <Link
                  className="stretched-link"
                  exact={exact}
                  key={index}
                  to={route}
                >
                  <span className="sr-only">{title + ' link'}</span>
                  <div
                    className="rounded-circle"
                    style={{
                      marginTop: '-40px',
                      padding: '10px',
                      backgroundColor: colors[index],
                      border: '4px solid white',
                      fillOpacity: '0.5',
                    }}
                  >
                    <img alt="icon" src={image} height="55" width="55" />
                  </div>
                </Link>
                <Card.Body>
                  <Card.Title className="text-dark">
                    <h2 style={{ fontSize: '1.75rem' }}>
                      <b>{cardTitle}</b>
                    </h2>
                  </Card.Title>
                  <Card.Text className="text-dark">
                    <small>{cardText}</small>
                  </Card.Text>
                </Card.Body>
              </Card>
              <div className="d-lg-none w-100"></div>
            </>
          )
        )}
      </CardDeck>
    );
  }

  return (
    <>
      <div className="banner-container text-center d-none d-md-block">
        <div className="image-blurred-edge">
          <img
            src="assets/images/banner-sample.jpg"
            alt="mSigPortal banner"
            style={{
              width: '100%',
              height: '250px',
            }}
          ></img>
        </div>
        <div className="banner-overlay-text row justify-content-center text-center text-light w-75">
          <div className="col-12">
            <h1 className="text-light">
              <b>mSigPortal</b>
            </h1>
          </div>
          <div
            className="col-6 w-50 my-3 align-self-center"
            style={{ borderTop: '3px solid white' }}
          ></div>
          <div
            className="col-12 text-center mt-2 font-weight-bold"
            style={{ width: '100%', fontSize: '18pt' }}
          >
            App Description
          </div>
          <div
            className="col-12 text-center mt-5"
            style={{ width: '100%', fontSize: '14pt' }}
          ></div>
        </div>
      </div>

      {/* mobile */}
      <div className="text-center mt-2 d-md-none">
        <h1 className="text-dark">
          <b>mSigPortal</b>
        </h1>
        <hr className="w-75"></hr>
        <div className="px-3 text-center">App Description</div>
      </div>

      <div
        className="container align-middle text-center"
        style={{ marginTop: '70px' }}
      >
        {cardRow(links.slice(0, 2), colors.slice(0, 2))}
        {cardRow(links.slice(2, 4), colors.slice(2, 4))}
      </div>
    </>
  );
}
