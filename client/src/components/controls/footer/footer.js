import React from 'react';
import { useLocation } from 'react-router-dom';

export function Footer() {
  const location = useLocation();
  return (
    <footer className="py-2 bg-gray">
      <div className={location.pathname === '/' ? '' : 'pl-2 pr-4'}>
        <div className="mt-4 text-light text-center">
          <div
            className={
              location.pathname === '/' ? 'container' : 'container-fluid '
            }
          >
            <div className="footer-nav text-left">
              <div className="mb-4 row justify-content-between">
                <div className="col-md-auto footer-nav-col">
                  <div className="d-none d-lg-block footer-header text-left mb-2">
                    <h5>
                      <b>Division of Cancer Epidemiology and Genetics</b>{' '}
                    </h5>
                    <p>
                      <b>at the National Cancer Institute</b>
                    </p>
                  </div>
                  <div className="d-lg-none footer-header-mobile text-left mb-2">
                    <h5>Division of Cancer Epidemiology and Genetics</h5>
                    <span>at the National Cancer Institute</span>
                  </div>
                </div>
              </div>
              <div className="row justify-content-between">
                <div className="col-md-auto footer-nav-col">
                  <h6 className="mb-2">CONTACT INFORMATION</h6>
                  <div className="my-2">
                    <a
                      className="footer-link text-light"
                      href="mailto:NCImSigPortalWebAdmin@mail.nih.gov"
                      target="_blank"
                    >
                      Contact Us
                    </a>
                  </div>
                </div>
                <div className="col-md-auto footer-nav-col">
                  <h6 className="mb-2">MORE INFORMATION</h6>
                  <div className="my-2">
                    <a
                      className="footer-link text-light"
                      href="https://dceg.cancer.gov/"
                      target="_blank"
                    >
                      DCEG
                    </a>
                  </div>
                  <div className="my-2">
                    <a
                      className="footer-link text-light"
                      href="https://prevention.cancer.gov/"
                      target="_blank"
                    >
                      DCP
                    </a>
                  </div>
                  {/* <!-- <div className="my-2">
                <a
                  className="footer-link text-light"
                  href="https://dceg.cancer.gov/research/who-we-study/cohorts/prostate-lung-colon-ovary-prospective-study"
                  target="_blank"
                  >PLCO Website</a
                >
              </div> --> */}
                </div>
                <div className="col-md-auto footer-nav-col">
                  <h6 className="mb-2">POLICIES</h6>
                  <div className="my-2">
                    <a
                      className="footer-link text-light"
                      href="https://www.cancer.gov/policies/accessibility"
                      target="_blank"
                    >
                      Accessibility
                    </a>
                  </div>
                  <div className="my-2">
                    <a
                      className="footer-link text-light"
                      href="https://www.cancer.gov/policies/disclaimer"
                      target="_blank"
                    >
                      Disclaimer
                    </a>
                  </div>
                  <div className="my-2">
                    <a
                      className="footer-link text-light"
                      href="https://www.cancer.gov/policies/foia"
                      target="_blank"
                    >
                      FOIA
                    </a>
                  </div>
                  <div className="my-2">
                    <a
                      className="footer-link text-light"
                      href="https://www.cancer.gov/policies/comments"
                      target="_blank"
                    >
                      Comment Policy
                    </a>
                  </div>
                  <div className="my-2">
                    <a
                      className="footer-link text-light"
                      href="https://www.hhs.gov/vulnerability-disclosure-policy/index.html"
                      target="_blank"
                    >
                      HHS Vulnerability Disclosure
                    </a>
                  </div>
                </div>
              </div>
            </div>

            <div className="footer-agencies my-4">
              <div className="row justify-content-md-center">
                <div className="col-sm-auto">
                  <a
                    className="footer-link text-light"
                    href="http://www.hhs.gov/"
                    target="_blank"
                  >
                    U.S. Department of Health and Human Services
                  </a>
                </div>
                <div className="col-sm-auto">
                  <a
                    className="footer-link text-light"
                    href="http://www.nih.gov"
                    target="_blank"
                  >
                    National Institutes of Health
                  </a>
                </div>
                <div className="col-sm-auto">
                  <a
                    className="footer-link text-light"
                    href="https://www.cancer.gov/"
                    target="_blank"
                  >
                    National Cancer Institute
                  </a>
                </div>
                <div className="col-sm-auto">
                  <a
                    className="footer-link text-light"
                    href="http://usa.gov"
                    target="_blank"
                  >
                    USA.gov
                  </a>
                </div>
              </div>
            </div>

            <div className="footer-tagline">
              <h6>
                NIH ... Turning Discovery Into Health
                <sup>®</sup>
              </h6>
            </div>
          </div>
        </div>
      </div>
    </footer>
  );
}
