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
              <div
                className="row row-cols-1 row-cols-md-3 g-4"
                key={title}
                style={{ marginRight: '1%' }}
              >
                <div class="col mb-4">
                  <Card key={title} id={title} className="card h-100 p-3">
                    {/* <Link
                      className="stretched-link"
                      exact={exact}
                      key={index}
                      to={route}
                    >
                      <span className="sr-only">{title + ' link'}</span>
                    </Link> */}
                    <img
                      alt={cardTitle}
                      src={image}
                      className="card-img-top w-40 h-40 ml-3"
                    />
                    <Card.Body>
                      <h4 class="card-title text-catalog">{cardTitle}</h4>
                      <p class="card-text">
                        {description}
                        <br />
                        <a href="/#/about" class="link-primary">
                          Read More &rarr;
                        </a>
                      </p>
                    </Card.Body>
                    <div class="p-3">
                      <a
                        href={route}
                        class="btn btn-lg btn-catalog btn-border-radius-25"
                      >
                        Go to {title} &gt;
                      </a>
                    </div>
                  </Card>
                </div>
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
        <div className="background-img text-center">
          <img
            src="assets/images/Hero_Image.png"
            alt="mSigPortal banner"
            style={{
              width: '100%',
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
                <h5>Learn more about mSigportal &gt;</h5>
              </div>
            </div>
          </div>
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

      <div className="p-3 home-grid-banner">
        {cardRow(links.slice(0, 1))}
        {cardRow(links.slice(1, 2))}
        {cardRow(links.slice(2, 3))}
        {cardRow(links.slice(3, 4))}
      </div>
    </>
  );
}
