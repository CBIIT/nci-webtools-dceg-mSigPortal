import React, { useState, useEffect } from 'react';
import { Link } from 'react-router-dom';
import Card from 'react-bootstrap/Card';
import { CardDeck, Button } from 'react-bootstrap';
import './home.scss';
import { ConsoleTransportOptions } from 'winston/lib/winston/transports';

export default function Home({ links }) {
  const colors = ['#fc8701', '#2c71dd', '#689f39', '#84368d'];

  function cardRow(links, color) {
    return (
      <CardDeck>
        {links.map(
          (
            { exact, route, action, title, cardTitle, cardText, description, image },
            index
          ) => (

              <div className="d-flex bd-highlight w-100" key={title} style={{ marginRight: '1%' }}>

                <Card
                  key={title}
                  id={title}
                  className="mb-5 align-self-center p-2 bd-highlight"
                  style={{
                    minWidth: '45%',
                    justifyContent: 'center',
                    position: 'relative',
                    border: '1px solid #DADBE6',
                    backgroundColor: 'white'
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

                  </Link>

                  <Card.Body>
                    <div className="d-flex flex-row align-items-center justify-content-left">
                      <div
                        className="rounded-circle"
                        style={{
                          //marginTop: '-15px',
                          marginLeft: '-10px',
                          marginRight: '10px',
                          width: '120px',
                          height: '120px',
                          padding: '10px',
                          backgroundColor: color,
                          border: '4px solid white',
                          fillOpacity: '0.5'
                        }}
                      >
                        <img alt="icon" src={image} height="105" width="105"
                          style={{
                            marginTop: '-13px',
                            marginLeft: '-7px'
                          }}
                        />

                      </div>
                      <Card.Title className="text-dark">
                        <h2 style={{ fontSize: '1.75rem', marginBottom:'-5px' }}>
                          <b>{cardTitle}</b>
                        </h2>
                      </Card.Title>
                    </div>
                    <div className="text-center mt-2 d-md-none">
                      {description}
                    </div>
                  </Card.Body>
                </Card>
                <div className="description d-none d-md-block">
                  {description}
                </div>

              </div>
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
            src="assets/images/msigportal-banner.png"
            alt="mSigPortal banner"
            style={{
              width: '100%',
              height: '250px',
              filter: 'contrast(2)'
            }}
          ></img>
        </div>
        <div className="banner-overlay-text row justify-content-center text-center text-light w-75">
          <div className="col-12">
            <img src="assets/images/logo-horizontal.png"
              style={{
                width: '325px',
                height: '50px'
              }}
            ></img>
          </div>
          <div
            className="col-6 w-50 my-3 align-self-center"
            style={{ borderTop: '3px solid', color: 'rgb(200,37,6)' }}
          ></div>
          <div
            className="col-12 text-center mt-2 font-weight-bold"
            style={{ width: '100%', fontSize: '18pt', color: 'black', fontStyle: 'italic' }}
          >
            Integrative mutational signature portal for cancer genomic study
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
        <div className="px-3 text-center">
          <b>Integrative mutational signature portal for cancer genomic study</b>
        </div>
      </div>

      <div
        className="align-middle"
        style={{
          marginLeft: '6%',
          marginRight: '6%',
          textAlign: 'left'
        }}
      >
        Mutational signatures are characteristic combination of mutation types arising from specific mutagenesis processes such as DNA replication infidelity, exogenous or endogenous exposures, defective DNA repair and DNA enzymatic editing. Analysis of mutational signatures is becoming routine in cancer genomics and provide a novel opportunity for biomarker discovery, tumor diagnostics and treatment guidance. mSigPortal will provide the state-of-the art methods and platform to explore, visualize and analyze mutational signatures, which will greatly facilitate broad investigation of mutational signatures to elucidate different mutagenesis processes involved in tumorigenesis.  Currently, mSigPortal includes the following modules:
      </div>

      <div
        className="container align-middle text-center"
        style={{ marginTop: '70px' }}
      >

        {cardRow(links.slice(0, 1), colors[0])}
        {cardRow(links.slice(1, 2), colors[1])}
        {cardRow(links.slice(2, 3), colors[2])}
        {cardRow(links.slice(3, 4), colors[3])}
      </div>
    </>
  );
}
