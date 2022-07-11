import React from 'react';
import { Link } from 'react-router-dom';
import { Button, OverlayTrigger, Popover } from 'react-bootstrap';
import { Animated } from 'react-animated-css';
import parse from 'html-react-parser';
import './home.scss';

export default function Home({ links }) {
  return (
    <>
      <div
        style={{
          backgroundImage: `url("assets/images/placeholder-home-banner-bg.png")`,
          backgroundRepeat: 'no-repeat',
          backgroundPosition: 'right',
        }}
      >
        <div className="container py-5">
          <div className="row my-5">
            <div className="col-12">
              <span
                style={{
                  fontSize: '4em',
                  color: '#D62D4C',
                  letterSpacing: '3px',
                }}
              >
                mSigPortal
              </span>
            </div>
            <div className="col-12 mt-1 mb-3">
              <span
                style={{
                  fontSize: '1.2em',
                  color: '#AA4A6D',
                  lineHeight: '1px',
                  letterSpacing: '1px',
                }}
              >
                Integrative Mutational Signature Portal
              </span>
              <br />
              <span
                style={{
                  fontSize: '1.2em',
                  color: '#2F54A5',
                  lineHeight: '1px',
                }}
              >
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
                    fontWeight: '500',
                  }}
                >
                  Explore the Catalog
                </Link>
              </Button>
            </div>
          </div>
        </div>
      </div>
      <div style={{ backgroundColor: '#283E5A', position: 'relative' }}>
        <div
          className="p-5 d-none d-md-block"
          style={{
            position: 'absolute',
            // right: '0px',
            height: '100%',
            overflow: 'hidden',
            width: '100%',
          }}
        >
          <Animated
            animationIn="slideInRight"
            animationInDuration={2000}
            isVisible={true}
            style={{ height: '100%' }}
          >
            <img
              alt={'DNA accent banner'}
              src={'assets/images/placeholder-home-dna.png'}
              height="100%"
              style={{ float: 'right' }}
            />
          </Animated>
        </div>
        <div className="container py-5">
          <div className="row mt-3 mb-5">
            <div className="col-12">
              <div className="home-text-border mb-4"></div>
            </div>
            <div className="col-12 text-light">
              <p
                style={{
                  fontWeight: '200',
                  lineHeight: '1.8em',
                }}
              >
                Mutational signatures are characteristic combinations of types
                arising from <br className="d-none d-md-block" />
                specific mutagenesis processes such as DNA replication
                infidelity, defective <br className="d-none d-md-block" />
                DNA repair, DNA enzymatic editing and exogenous exposures.
                Analysis of <br className="d-none d-md-block" />
                mutational signatures is becoming routine in cancer genomics,
                providing a <br className="d-none d-md-block" />
                novel opportunity for biomarker discovery, tumor diagnostics,
                and treatment <br className="d-none d-md-block" />
                guidance. As the number of mutational signatures associated with
                known <br className="d-none d-md-block" />
                etiologies has increased from many different cancer genomic
                studies, there is <br className="d-none d-md-block" />a critical
                need for curated census as well as data sharing of mutational{' '}
                <br className="d-none d-md-block" />
                signatures for public research. mSigPortal provides a platform
                that enables <br className="d-none d-md-block" />
                users to explore, visualize, and analyze mutational signatures
                and relevant <br className="d-none d-md-block" />
                signature data (such as mutational profile, proposed etiology,
                tissue <br className="d-none d-md-block" />
                specificity, activity, and association) in cancer genomic
                studies from scientific <br className="d-none d-md-block" />
                literature or user input.{' '}
              </p>
            </div>
          </div>
        </div>
      </div>
      <div className="bg-white py-5">
        <div className="mx-5 my-4 px-2">
          <div className="row justify-content-center">
            {links
              .slice(0, 4)
              .map(
                (
                  { cardTitle, image, exact, route, title, description },
                  idx
                ) => (
                  <OverlayTrigger
                    placement="top"
                    // delay={{ hide: 1000 }}
                    overlay={(props) => (
                      <Popover {...props} className="p-3">
                        <p>{parse(description)}</p>
                      </Popover>
                    )}
                  >
                    <div
                      className="col-auto d-flex align-items-center mx-1 my-1 home-nav-card"
                      key={`home-nav-card-${idx}`}
                      style={{
                        backgroundColor: '#EFEFEF',
                        width: '300px',
                        height: '300px',
                      }}
                    >
                      <Link className="stretched-link" exact={exact} to={route}>
                        <span className="sr-only">{title + ' link'}</span>
                      </Link>
                      <div className="text-center w-100">
                        <div className="">
                          <img
                            alt={title}
                            src={image}
                            height="150"
                            width="150"
                          />
                        </div>
                        <div className="mt-3">
                          {cardTitle.split(' ').map((titleItem, idx) => (
                            <div
                              key={`home-nav-card-title-item-${titleItem}-${idx}`}
                            >
                              <span
                                style={{
                                  fontWeight: '800',
                                  color: '#366D83',
                                  letterSpacing: '1px',
                                }}
                              >
                                {titleItem}
                              </span>
                            </div>
                          ))}
                        </div>
                      </div>
                    </div>
                  </OverlayTrigger>
                )
              )}
          </div>
        </div>
      </div>
    </>
  );
}
