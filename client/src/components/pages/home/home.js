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
        <div className="">
          <img
            src="/assets/images/Hero_Image.png"
            alt="mSigPortal banner"
            style={{
              width: '100%',
              //height: '250px',
              //filter: 'contrast(2)',
            }}
          ></img>
          <div class="row">
            <div class="homepage-title-left text-left">
              <h1>mSigPortal</h1>
              <div class="text-primary-purple msigportal-title">
                Integrative Mutational Signature Portal for Cancer Genomics
                Study
              </div>
            </div>
            <div class="homepage-title-right">
              <div class=" btn btn-gradient btn-1 p-3">
                Learn more about mSigportal &gt;
              </div>
            </div>
          </div>
        </div>
      </div>
      <div className=" p-3">
        <div class="row row-cols-1 row-cols-md-3 g-4">
          <div class="col mb-4">
            <div class="card h-100 p-3">
              <img
                src="/assets/icons/Catalog-Icon.svg"
                class="card-img-top w-50 h-50 ml-3"
                alt="Catalog Icon"
              />
              <div class="card-body">
                <h4 class="card-title text-catalog">Signature Catalog</h4>
                <p class="card-text">
                  All existing human and mouse signatures based on different
                  genome builds and algorithm versions <br />
                  <a href="/#/about">Read More</a>
                </p>
              </div>
              <div class="p-3">
                <a
                  href="/mutational-signatures/#/catalog"
                  class="btn btn-catalog btn-border-radius-25"
                >
                  Go to catalog &gt;
                </a>
              </div>
            </div>
          </div>
          <div class="col mb-3">
            <div class="card h-100 p-3">
              <img
                src="/assets/icons/Visualization-Icon.svg"
                class="card-img-top w-50 h-50 ml-3"
                alt="Visualization Icon"
              />
              <div class="card-body">
                <h4 class="card-title text-visualization">
                  Signature Visualization
                </h4>
                <p class="card-text">
                  Allow identication of signature features at sample level and
                  sicovery of new signatures <br />
                  <a href="/mutational-signatures/#/about">Read More</a>
                </p>
              </div>
              <div class="p-3">
                <a
                  href="/mutational-signatures/#/visualization"
                  class="btn btn-visualization"
                >
                  Go to visualization &gt;
                </a>
              </div>
            </div>
          </div>
          <div class="col mb-3">
            <div class="card h-100 p-3">
              <img
                src="/assets/icons/Extraction-Icon.svg"
                class="card-img-top w-50 h-50 ml-3"
                alt="Extraction Icon"
              />
              <div class="card-body">
                <h4 class="card-title text-extraction">Signature Extraction</h4>
                <p class="card-text">
                  Extract and compare muational signatures using
                  state-of-the-art algorithms <br />
                  <a href="/mutational-signatures/#/about">Read More</a>
                </p>
              </div>
              <div class="p-3">
                <a
                  href="/mutational-signatures/#/extraction"
                  class="btn btn-extraction"
                >
                  Go to extraction &gt;
                </a>
              </div>
            </div>
          </div>
          <div class="col mb-3">
            <div class="card h-100 p-3">
              <img
                src="/assets/icons/Exploration-Icon.svg"
                class="card-img-top w-50 h-50 ml-3"
                alt="Exploration Icon"
              />
              <div class="card-body">
                <h4 class="card-title text-exploration">
                  Signature Exploration
                </h4>
                <p class="card-text">
                  Explore etiological factors associated with signature at
                  sample level <br />
                  <a href="/mutational-signatures/#/about">Read More</a>
                </p>
              </div>
              <div class="p-3">
                <a
                  href="/mutational-signatures/#/exploration"
                  class="btn btn-exploration"
                >
                  Go to exploration &gt;
                </a>
              </div>
            </div>
          </div>
          <div class="col mb-3">
            <div class="card h-100 p-3">
              <img
                src="/assets/icons/Association-Icon.svg"
                class="card-img-top  w-50 h-50 ml-3"
                alt="Association Icon"
              />
              <div class="card-body">
                <h4 class="card-title text-association">
                  Signature Association
                </h4>
                <p class="card-text">
                  Analyze signature association with other genomic features and
                  clincial data <br />
                  <a href="/mutational-signatures/#/about">Read More</a>
                </p>
              </div>
              <div class="p-3">
                <a
                  href="/mutational-signatures/#/association"
                  class="btn btn-association"
                >
                  Go to Association &gt;
                </a>
              </div>
            </div>
          </div>
          <div class="col mb-3">
            <div class="card h-100 p-3">
              <img
                src="/assets/icons/API-Icon.svg"
                class="card-img-top w-50 h-50 ml-3"
                alt="API Icon"
              />
              <div class="card-body">
                <h4 class="card-title text-api">Signature API Access</h4>
                <p class="card-text">
                  Analyze signature association with other genomic features and
                  clincial data <br />
                  <a href="/mutational-signatures/#/about">Read More</a>
                </p>
              </div>
              <div class="p-3">
                <a
                  href="/mutational-signatures/#/apiaccess"
                  class="btn btn-api"
                >
                  Go to API Access &gt;
                </a>
              </div>
            </div>
          </div>
        </div>
      </div>
    </>
  );
}
