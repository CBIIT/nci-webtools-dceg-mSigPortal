import React from 'react';
import { Link } from 'react-router-dom';
import Card from 'react-bootstrap/Card';
import { Button, CardDeck } from 'react-bootstrap';
import parse from 'html-react-parser';
import './home.scss';
export default function Home({ links }) {
  function cardRow(links) {
    return (
      <div>
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
              <div key={title}>
                <Card
                  key={title}
                  id={title}
                  className=""
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
                    <div className="">
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
                        <h2
                          style={{ fontSize: '1.75rem', marginBottom: '-5px' }}
                        >
                          <b>{cardTitle}</b>
                        </h2>
                      </Card.Title>
                      <div className="description d-none d-md-block">
                        <div>{parse(description)}</div>
                      </div>
                    </div>
                  </Card.Body>
                </Card>
              </div>
            )
          )}
        </CardDeck>
      </div>
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
        className="container mb-4"
        // style={{ marginTop: '70px' }}
      >
        <div class="row row-cols-1 row-cols-md-3 g-4">
          <div class="col mb-3">
            <div class="card h-100">
              <img src="..." class="card-img-top" alt="..." />
              <div class="card-body">
                <h4 class="card-title text-catalog">Signature Catalog</h4>
                <p class="card-text">
                  All existing human and mouse signatures based on different
                  genome builds and algorithm versions <br />
                  <a href="/#/about">Read More</a>
                </p>
              </div>
              <div class="p-3">
                <a href="/#/catalog" class="btn btn-catalog">
                  Go to catalog
                </a>
              </div>
            </div>
          </div>
          <div class="col mb-3">
            <div class="card h-100">
              <img src="..." class="card-img-top" alt="..." />
              <div class="card-body">
                <h4 class="card-title text-visualization">
                  Signature Visualization
                </h4>
                <p class="card-text">
                  Allow identication of signature features at sample level and
                  sicovery of new signatures <br />
                  <a href="/#/about">Read More</a>
                </p>
              </div>
              <div class="p-3">
                <a href="/#/visualization" class="btn btn-visualization">
                  Go to visualization
                </a>
              </div>
            </div>
          </div>
          <div class="col mb-3">
            <div class="card h-100">
              <img src="..." class="card-img-top" alt="..." />
              <div class="card-body">
                <h4 class="card-title text-extraction">Signature Extraction</h4>
                <p class="card-text">
                  Extract and compare muational signatures using
                  state-of-the-art algorithms <br />
                  <a href="/#/about">Read More</a>
                </p>
              </div>
              <div class="p-3">
                <a href="/#/extraction" class="btn btn-extraction">
                  Go to extraction
                </a>
              </div>
            </div>
          </div>
          <div class="col mb-3">
            <div class="card h-100">
              <img src="..." class="card-img-top" alt="..." />
              <div class="card-body">
                <h4 class="card-title text-exploration">
                  Signature Exploration
                </h4>
                <p class="card-text">
                  Explore etiological factors associated with signature at
                  sample level <br />
                  <a href="/#/about">Read More</a>
                </p>
              </div>
              <div class="p-3">
                <a href="/#/exploration" class="btn btn-exploration">
                  Go to exploration
                </a>
              </div>
            </div>
          </div>
          <div class="col mb-3">
            <div class="card h-100">
              <img src="..." class="card-img-top" alt="..." />
              <div class="card-body">
                <h4 class="card-title text-association">
                  Signature Association
                </h4>
                <p class="card-text">
                  Analyze signature association with other genomic features and
                  clincial data <br />
                  <a href="/#/about">Read More</a>
                </p>
              </div>
              <div class="p-3">
                <a href="/#/association" class="btn btn-association">
                  Go to Association
                </a>
              </div>
            </div>
          </div>
          <div class="col mb-3">
            <div class="card h-100">
              <img src="..." class="card-img-top" alt="..." />
              <div class="card-body">
                <h4 class="card-title text-api">Signature API Access</h4>
                <p class="card-text">
                  Analyze signature association with other genomic features and
                  clincial data <br />
                  <a href="/#/about">Read More</a>
                </p>
              </div>
              <div class="p-3">
                <a href="/#/apiaccess" class="btn btn-api">
                  Go to API Access
                </a>
              </div>
            </div>
          </div>
        </div>
      </div>
    </>
  );
}
