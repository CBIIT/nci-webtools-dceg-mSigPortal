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
        <div className="col">
          <Card
            key={title}
            id={title}
            className="p-3"
            style={{
              borderRadius: '2em',
              display: 'inline-block',
            }}
          >
            <img
              alt={cardTitle}
              src={image}
              className="card-img-top w-25 h-25 ml-3"
            />
            <Card.Body>
              <h4 className={`card-title card-title-homepage text-${name}`}>
                {cardTitle}
              </h4>
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
              <div class="">
                <Link
                  className={`btn btn-2 btn-${name} btn-border-radius-15`}
                  exact={exact}
                  key={index}
                  to={route}
                >
                  <span class="">Go to {title} &gt;</span>
                </Link>
              </div>
            </Card.Body>
          </Card>
        </div>
      </div>
    );
  }

  function CardRow2(
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
      <div key={index} className="col-sm-12 col-md-6 col-lg-4 mb-3" id={title}>
        <div
          style={{
            borderRadius: '2em',
            display: 'inline-block',
            backgroundColor: 'white',
          }}
        >
          <div className="flex-fill">
            <div className="m-3">
              <img
                alt={cardTitle}
                src={image}
                className="card-img-top w-25 h-25"
              />
              <h4 className={`card-title card-title-homepage text-${name}`}>
                {cardTitle}
              </h4>
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
              <div class="">
                <Link
                  className={`btn btn-2 btn-${name} btn-border-radius-15`}
                  exact={exact}
                  key={index}
                  to={route}
                >
                  <span class="">Go to {title} &gt;</span>
                </Link>
              </div>
            </div>
          </div>
        </div>
      </div>
    );
  }

  function CardRow3(
    {
      exact,
      route,
      action,
      title,
      cardId,
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
      <div
        key={index}
        class="col-lg-4 col-md-6 col-sm-12 mb-3 d-flex align-items-stretch"
      >
        <div
          className="card"
          id={cardId}
          style={{
            borderRadius: '2em',
            display: 'inline-block',
            backgroundColor: 'white',
          }}
        >
          <img
            alt={cardTitle}
            src={image}
            className="card-img-top w-40 ml-3 mt-3"
          />
          <div className="card-body d-flex flex-column">
            <h5 className={`card-title card-title-homepage text-${name}`}>
              {cardTitle}
            </h5>
            <p class="card-text mb-4">
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
            <div class="pb-1">
              <Link
                className={`btn btn-2 btn-${name} btn-border-radius-15`}
                exact={exact}
                key={index}
                to={route}
              >
                <span
                  style={{
                    fontSize: '13px',
                  }}
                >
                  Go to {title} &gt;
                </span>
              </Link>
            </div>
          </div>
        </div>
      </div>
    );
  }

  return (
    <>
      <div className="banner-container">
        <div className="image-blurred-edge">
          <div className="background-img text-center">
            <img
              src="assets/images/Hero_Image.png"
              alt="mSigPortal banner"
              style={{
                width: '100%',
                height: '365px',
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
              <div className="homepage-title-right">
                <div className=" btn btn-gradient btn-1">
                  <div className="msigportal-home-title-right">
                    Learn more about mSigportal &gt;
                  </div>
                </div>
              </div>
            </div>
          </div>
        </div>

        <div
          className="card-grid"
          style={{
            zIndex: '1',
            position: 'relative',
          }}
        >
          {/* <div className="container">
            <div className="row">
              {links
                .filter((e) => e.showHomepage)
                .map((e, i) => CardRow3(e, i))}
            </div>
          </div> */}
          <div className="row">
            <div className="col-lg-2 col-md-1 col-sm-1 col-xs-1"></div>
            <div className="col-lg-8 col-md-10 col-sm-10 col-xs-10">
              <div className="row">
                {links
                  .filter((e) => e.showHomepage)
                  .map((e, i) => CardRow3(e, i))}
              </div>
            </div>
            <div className="col-lg-2 col-md-1 col-sm-1 col-xs-1"></div>
          </div>
        </div>
      </div>
    </>
  );
}
