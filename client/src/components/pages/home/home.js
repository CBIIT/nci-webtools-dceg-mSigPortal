import React, { useState, useEffect } from 'react';
import { Link } from 'react-router-dom';
import Card from 'react-bootstrap/Card';
import { CardDeck, Button } from 'react-bootstrap';
import './home.scss';

export default function Home({ links }) {
  const colors = ['#ffbe63', '#799eff', '#a4ce76', '#ab68b9'];

  function cardRow(links) {
    return (
      <CardDeck>
        {links.map(
          (
            { exact, route, action, title, cardTitle, cardText, description, image },
            index
          ) => (

            <div className="d-flex bd-highlight w-100" key={title}>

              <Card
                key={title}
                id={title}
                className="mb-5 align-self-center p-2 bd-highlight"
                style={{
                  minWidth: '45%',
                  justifyContent: 'center',
                  alignItems: 'center',
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
                  <img alt="icon" src={image} height="80" width="80" 
                    style={{
                      marginTop: '-40px',
                      borderRadius: '80'
                    }} />
                  {/*
                  <span className="sr-only">{title + ' link'}</span>
                  <div
                    className="rounded-circle"
                    style={{
                      marginTop: '-40px',
                      padding: '10px',
                      backgroundColor: color,
                      border: '4px solid white',
                      fillOpacity: '0.5'
                    }}
                  >
                    
                  </div>*/}
                </Link>

                <Card.Body>
                  <Card.Title className="text-dark">
                    <h2 style={{ fontSize: '1.75rem' }}>
                      <b>{cardTitle}</b>
                    </h2>
                   </Card.Title>
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
            <h1 
              style={{color:'rgb(200,37,6)', font: '40px Phosphate Inline, sans-serif', fontWeight:'300'}}
            >
              <b>mSigPortal</b>
            </h1>
          </div>
          <div
            className="col-6 w-50 my-3 align-self-center"
            style={{ borderTop: '3px solid white' }}
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
        className="container align-middle "
        style={{ 
          marginTop: '50px', 
          textAlign: 'left'
        }}
      >
        Mutational signatures are characteristic combination of mutation types arising from specific mutagenesis processes such as DNA replication infidelity, exogenous or endogenous exposures, defective DNA repair and DNA enzymatic editing. Analysis of mutational signatures is becoming routine in cancer genomics and provide a novel opportunity for biomarker discovery, tumor diagnostics and treatment guidance. mSigPortal will provide the state-of-the art methods and platform to explore, visualize and analyze mutational signatures, which will greatly facilitate broad investigation of mutational signatures to elucidate different mutagenesis processes involved in tumorigenesis.  Currently, mSigPortal includes the following modules:
      </div>
        
      <div
        className="container align-middle text-center"
        style={{ marginTop: '70px' }}
      >
        
        {cardRow(links.slice(0, 1))}
        {cardRow(links.slice(1, 2))}
        {cardRow(links.slice(2, 3))}
        {cardRow(links.slice(3, 4))}
      </div>
    </>
  );
}
