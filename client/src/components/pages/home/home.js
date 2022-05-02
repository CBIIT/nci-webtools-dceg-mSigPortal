import React from 'react';
import { Link } from 'react-router-dom';
import Card from 'react-bootstrap/Card';
import { CardDeck, Button } from 'react-bootstrap';
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
            },
            index
          ) => (
            <div
              className="d-flex bd-highlight w-100 mb-5"
              key={title}
              style={{ marginRight: '1%' }}
            >
              <Card
                key={title}
                id={title}
                className="align-self-center p-2 bd-highlight"
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
                <p>{parse(description)}</p>
              </div>
            </div>
          )
        )}
      </CardDeck>
    );
  }

  return (
    <>
      <div style={{
        backgroundImage: `url("assets/images/placeholder-home-banner-bg.png")`,
        backgroundRepeat: 'no-repeat',
        backgroundPosition: 'right'
      }}>
        <div className="container py-5">
          <div className="row my-5">
            <div className="col-12">
              <span style={{
                fontSize: '4em',
                color: '#D62D4C',
                letterSpacing: '3px'
                }}>
                mSIGPORTAL
              </span>
            </div>
            <div className="col-12 mt-1 mb-3">
              <span style={{
                fontSize: '1.2em',
                color: '#AA4A6D',
                lineHeight: '1px',
                letterSpacing: '1px'
              }}>
                Integrative Mutational Signature Portal
              </span>
              <br />
              <span style={{
                fontSize: '1.2em',
                color: '#2F54A5',
                lineHeight: '1px'
              }}>
                for Cancer Genomic Studies
              </span>
            </div>
            <div className="col-12 my-2">
              <p>
                Facilitating broad investigation of mutational 
                <br />
                signatures to elucidate different mutagenesis 
                <br />
                processes involved in tumorigenesis.
              </p>
            </div>
            <div className="col-12 my-2">
              <Button
                variant={'outline-dark'}
                className="bg-white home-banner-nav-button"
                onClick={(e) => e.preventDefault()}
              >
                <Link
                  className="stretched-link m-2"
                  to={'/catalog'}
                  style={{
                    textDecoration: 'none',
                    color: '#D62D4C',
                    fontWeight: '500'
                  }}>
                  Explore the Catalog
                </Link>
              </Button>
            </div>
          </div>
        </div>
      </div>
      <div style={{backgroundColor: '#283E5A'}}>
        <div className="container py-5">
          <div className="row mt-3 mb-5">
            <div className="col-12">
              <div className='home-text-border mb-4'></div>
            </div>
            <div className="col-12 text-light">
              <p style={{
                fontWeight: '200',
                lineHeight: '1.8em'
              }}>
                Mutational signatures are characteristic combinations of types arising from{' '}
                <br className="d-none d-md-block" />
                specific mutagenesis processes such as DNA replication infidelity, defective{' '}
                <br className="d-none d-md-block" />
                DNA repair, DNA enzymatic editing and exogenous exposures. Analysis of{' '}
                <br className="d-none d-md-block" />
                mutational signatures is becoming routine in cancer genomics, providing a{' '} 
                <br className="d-none d-md-block" />
                novel opportunity for biomarker discovery, tumor diagnostics, and treatment{' '} 
                <br className="d-none d-md-block" />
                guidance. As the number of mutational signatures associated with known{' '} 
                <br className="d-none d-md-block" />
                etiologies has increased from many different cancer genomic studies, there is{' '} 
                <br className="d-none d-md-block" />
                a critical need for curated census as well as data sharing of mutational{' '} 
                <br className="d-none d-md-block" />
                signatures for public research. mSigPortal provides a platform that enables{' '} 
                <br className="d-none d-md-block" />
                users to explore, visualize, and analyze mutational signatures and relevant{' '} 
                <br className="d-none d-md-block" />
                signature data (such as mutational profile, proposed etiology, tissue{' '} 
                <br className="d-none d-md-block" />
                specificity, activity, and association) in cancer genomic studies from scientific{' '}
                <br className="d-none d-md-block" />
                literature or user input.{' '}
              </p>
            </div>
          </div>
        </div>
      </div>
      {/* <div className="banner-container text-center d-none d-md-block">
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
      </div> */}

      {/* mobile */}
      {/* <div className="text-center mt-2 d-md-none">
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
      </div> */}
    </>
  );
}
