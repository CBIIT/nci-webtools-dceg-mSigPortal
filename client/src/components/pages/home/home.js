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

  return (
    <>
      <div className="banner-container ">
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
        {/* <div className="home-grid-banner">
          <div className="home-grid-banner-container">
            <div className="row row-cols-1 row-cols-md-3 g-4 m-3">
              {links.filter((e) => e.showHomepage).map((e, i) => CardRow(e, i))}
            </div>
          </div>
        </div> */}

        <div className="home-grid-banner">
          <div className="home-grid-banner-container">
            <div className="card-columns">
              {links.filter((e) => e.showHomepage).map((e, i) => CardRow(e, i))}
            </div>
          </div>
        </div>
      </div>
    </>
  );
}
