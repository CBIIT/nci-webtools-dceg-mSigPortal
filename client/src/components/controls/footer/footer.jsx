import React from 'react';
import { useLocation } from 'react-router-dom';

function parseVersionAndDate(versionString) {
  if (!versionString)
    return { version: 'N/A', date: new Date().toISOString().split('T')[0] };
  const versionMatch = versionString.match(/(\d+\.\d+\.\d+)(_dev)?/);
  const version = versionMatch
    ? versionMatch[1] + (versionMatch[2] || '')
    : 'N/A';

  // Extract 8-digit date if present
  const dateMatch = versionString.match(/(\d{8})/)?.[1];
  const date = dateMatch
    ? `${dateMatch.slice(0, 4)}-${dateMatch.slice(4, 6)}-${dateMatch.slice(
        6,
        8
      )}`
    : new Date().toISOString().split('T')[0];

  return { version, date };
}

export function Footer() {
  const location = useLocation();
  const { version, date } = parseVersionAndDate(
    import.meta.env.VITE_APP_VERSION
  );

  return (
    <div className="py-2 bg-gray">
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
                  <div className="d-none d-lg-block footer-header text-left ">
                    <h2 className="h5">
                      <b>Division of Cancer Epidemiology and Genetics</b>{' '}
                    </h2>
                    <div className="h6">
                      at the
                      <a className="text-white" href="https://www.cancer.gov/">
                        {' '}
                        National Institutes of Health{' '}
                      </a>
                    </div>
                  </div>
                  <div className="d-lg-none footer-header-mobile text-left mb-2">
                    <h5>Division of Cancer Epidemiology and Genetics</h5>
                    <span>at the National Cancer Institute</span>
                  </div>
                </div>
              </div>
              <div className="row justify-content-between">
                <div className="col-md-auto footer-nav-col">
                  <h2 className="h6 mb-1">CONTACT INFORMATION</h2>
                  <div className="my-0">
                    <a
                      className="footer-link text-light"
                      href="mailto:NCImSigPortalWebAdmin@mail.nih.gov"
                      target="_blank"
                    >
                      Contact Us
                    </a>
                  </div>
                  <div className="my-0">
                    <span className="text-light">
                      <b>Version: </b>
                      {version}
                    </span>
                  </div>
                  <div className="my-0">
                    <span className="text-light">
                      <b>Last Update: </b>
                      {date}
                    </span>
                  </div>
                </div>
                <div className="col-md-auto footer-nav-col">
                  <h2 className="h6 mb-1">MORE INFORMATION</h2>
                  <div className="my-0">
                    <a
                      className="footer-link text-light"
                      href="https://dceg.cancer.gov/"
                      target="_blank"
                    >
                      DCEG
                    </a>
                  </div>
                  {/* <div className="">
                    <a
                      className="footer-link text-light"
                      href="https://prevention.cancer.gov/"
                      target="_blank"
                    >
                      DCP
                    </a>
                  </div> */}
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
                  <h2 className="h6 mb-1">POLICIES</h2>
                  <div className="my-0">
                    <a
                      className="footer-link text-light"
                      href="https://www.cancer.gov/policies/accessibility"
                      target="_blank"
                    >
                      Accessibility
                    </a>
                  </div>
                  <div className="my-0">
                    <a
                      className="footer-link text-light"
                      href="https://www.cancer.gov/policies/disclaimer"
                      target="_blank"
                    >
                      Disclaimer
                    </a>
                  </div>
                  <div className="my-0">
                    <a
                      className="footer-link text-light"
                      href="https://www.cancer.gov/policies/foia"
                      target="_blank"
                    >
                      FOIA
                    </a>
                  </div>
                  <div className="my-0">
                    <a
                      className="footer-link text-light"
                      href="https://www.cancer.gov/policies/comments"
                      target="_blank"
                    >
                      Comment Policy
                    </a>
                  </div>
                  <div className="my-0">
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

            <div className="footer-agencies my-3">
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
              {/* <h6> */}
              NIH ... Turning Discovery Into Health
              <sup>®</sup>
              {/* </h6> */}
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
