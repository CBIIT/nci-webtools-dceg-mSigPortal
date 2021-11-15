import React from 'react';
import { Link } from 'react-router-dom';
import Card from 'react-bootstrap/Card';
import { CardDeck } from 'react-bootstrap';
import parse from 'html-react-parser';
import './home.scss';

export default function Home({ links }) {
  function cardRow(links) {
    return (
      <CardDeck>
        {links.map(
          (
            {
              exact,
              route,
              action,
              title,
              cardTitle,
              cardText,
              description,
              image,
              color,
              examples,
            },
            index
          ) => (
            <div
              className="d-flex bd-highlight w-100"
              key={title}
              style={{ marginRight: '1%' }}
            >
              <Card
                key={title}
                id={title}
                className="mb-5 align-self-center p-2 bd-highlight"
                style={{
                  minWidth: '45%',
                  justifyContent: 'center',
                  position: 'relative',
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
                        fillOpacity: '0.5',
                      }}
                    >
                      <img
                        alt={cardTitle}
                        src={image}
                        height="105"
                        width="105"
                        style={{
                          marginTop: '-13px',
                          marginLeft: '-7px',
                        }}
                      />
                    </div>
                    <Card.Title className="text-dark">
                      <h2 style={{ fontSize: '1.75rem', marginBottom: '-5px' }}>
                        <b>{cardTitle}</b>
                      </h2>
                    </Card.Title>
                  </div>
                </Card.Body>
              </Card>
              <div className="description d-none d-md-block">
                <div>{parse(description)}</div>
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
              filter: 'contrast(2)',
            }}
          ></img>
        </div>
        <div className="banner-overlay-text row justify-content-center text-center text-light w-75">
          <div className="col-12">
            <img
              src="assets/images/logo-horizontal.png"
              alt="mSigPortal title"
              style={{
                width: '325px',
                height: '50px',
              }}
            ></img>
          </div>
          <div
            className="col-6 w-50 my-3 align-self-center"
            style={{ borderTop: '3px solid', color: 'rgb(200,37,6)' }}
          ></div>
          <div
            className="col-12 text-center mt-2 font-weight-bold"
            style={{
              width: '100%',
              fontSize: '18pt',
              color: 'black',
              fontStyle: 'italic',
            }}
          >
            Integrative Mutational Signature Portal for Cancer Genomic Studies
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
          <b>
            Integrative mutational signature portal for cancer genomic studies
          </b>
        </div>
      </div>

      <div
        className="container align-middle"
        style={{
          textAlign: 'left',
        }}
      >
        Mutational signatures are characteristic combinations of mutation types
        arising from specific mutagenesis processes such as DNA replication
        infidelity, defective DNA repair, DNA enzymatic editing and exogenous
        exposures. Analysis of mutational signatures is becoming routine in
        cancer genomics, providing a novel opportunity for biomarker discovery,
        tumor diagnostics, and treatment guidance. As the number of mutational
        signatures associated with known etiologies has increased from many
        different cancer genomic studies, there is a critical need for curated
        census as well as data sharing of mutational signatures for public
        research. mSigPortal provides a platform that enables users to explore,
        visualize, and analyze mutational signatures and relevant signature data
        (such as mutational profile, proposed etiology, tissue specificity,
        activity, and association) in cancer genomic studies from scientific
        literature or user input. This portal will greatly facilitate broad
        investigation of mutational signatures to elucidate different
        mutagenesis processes involved in tumorigenesis. Currently, mSigPortal
        includes the following modules:
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
