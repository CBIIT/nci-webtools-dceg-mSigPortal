import React from 'react';
import { useLocation } from 'react-router-dom';
import HeaderSearch from './header-search';
export function Header() {
  const { pathname } = useLocation();

  return (
    <div className="bg-light">
      <include-html src="https://cbiit.github.io/nci-softwaresolutions-elements/banners/government-shutdown.html"></include-html>
      <a
        href="/#/#root"
        className="
            sr-only sr-only-focusable
            d-block
            text-white
            bg-primary-dark
            text-center
          "
      >
        Skip to Main Content
      </a>
      <div className="usa-banner" aria-label="Official government website">
        <div
          className={`d-none d-md-block ${
            pathname === '/' ? 'container p-0' : ''
          }`}
        >
          <div className="usa-banner__header">
            <div className="usa-banner__inner">
              <div
                className="usa-banner__header-text"
                style={{ whiteSpace: 'nowrap' }}
              >
                An official website of the United States government
              </div>
            </div>
          </div>
        </div>
      </div>
      <div>
        <div
          className={`d-none d-md-block ${
            pathname === '/' ? 'container p-0' : ''
          }`}
        >
          <div className="d-flex justify-content-between">
            <a href="https://dceg.cancer.gov/" target="_blank">
              <img
                src="https://analysistools.cancer.gov/common/images/mSigPortal_Desktop-Logo-COLOR.svg"
                height="100"
                alt="National Cancer Institute Logo"
              />
            </a>
            <div className="mr-4 mt-4">
              <HeaderSearch />
            </div>
          </div>
        </div>

        <div className="d-block d-sm-block d-md-none">
          <a href="https://dceg.cancer.gov/">
            <img
              src="https://analysistools.cancer.gov/common/images/mSigPortal_Desktop-Logo-COLOR.svg"
              height="80"
              width="100%"
              alt="National Cancer Institute Logo"
            />
          </a>
        </div>
      </div>
    </div>
  );
}
