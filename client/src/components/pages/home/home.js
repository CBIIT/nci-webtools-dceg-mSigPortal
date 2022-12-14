import React from 'react';
import { Link } from 'react-router-dom';
import Card from 'react-bootstrap/Card';
import { Button, CardDeck } from 'react-bootstrap';
import parse from 'html-react-parser';
import './home.scss';
export default function Home({ links }) {
  function CardRow(
    {
      exact,
      route,
      action,
      title,
      name,
      cardTitle,
      cardText,
      description,
      image,
      about,
      color,
      examples,
      buttonHomepage,
    },
    index
  ) {
    return (
      <div key={index}>
        <div class="col mb-3 ">
          <Card
            key={title}
            id={title}
            className="h-100 p-3"
            style={{
              borderRadius: '2em',
            }}
          >
            <img
              alt={cardTitle}
              src={image}
              className="card-img-top w-25 h-25 ml-3"
            />
            <Card.Body>
              <h4 className={`card-title text-${name}`}>{cardTitle}</h4>
              <p className="card-text">
                {description}
                <br />

                <Link
                  className="link-primary-underline"
                  exact={exact}
                  key={index}
                  to={about}
                >
                  <span>Read More &rarr;</span>
                </Link>
              </p>
            </Card.Body>
            <div class="p-3">
              <Link
                className={`btn btn-2 btn-lg btn-${name} btn-border-radius-25`}
                exact={exact}
                key={index}
                to={route}
              >
                <span class="">Go to {title} &gt;</span>
              </Link>
            </div>
          </Card>
        </div>
      </div>
    );
  }

  return (
    <>
      <div className="banner-container d-none d-md-block">
        <div className="background-img text-center">
          <img
            src="assets/images/Hero_Image.png"
            alt="mSigPortal banner"
            style={{
              width: '100%',
            }}
          ></img>
          <div className="row">
            <div className="homepage-title-left text-left">
              <h1 className="msigportal-home-title">mSigPortal</h1>
              <div class="text-primary-purple msigportal-title">
                Integrative Mutational Signature Portal for Cancer Genomics
                Study
              </div>
            </div>
            <div class="homepage-title-right">
              <div class=" btn btn-gradient btn-1 p-3">
                <div>Learn more about mSigportal &gt;</div>
              </div>
            </div>
          </div>
        </div>
        <div className="p-3 home-grid-banner">
          <div class="row row-cols-1 row-cols-md-3 g-4">
            {links.filter((e) => e.showHomepage).map((e, i) => CardRow(e, i))}
          </div>
        </div>
      </div>
    </>
  );
}
