import React from 'react';

export function Header() {
  return (
    <header className="bg-light">
      <div>
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
        <div className="ml-4 mr-5">
          <div className="p-2 d-none d-sm-block">
            <a href="https://dceg.cancer.gov/" target="_blank">
              <img
                src="https://analysistools.cancer.gov/common/images/DCEG-logo.svg"
                height="100"
                alt="National Cancer Institute Logo"
              />
            </a>
            {/* <a href="/#/">
              <img
                className="d-none d-md-block"
                src="assets/images/msigportal-logo.png"
                alt="mSigPortal Logo"
                height="80"
                style={{float: 'right', marginTop: '15px'}}
              />
            </a> */}
          </div>
          <div className="p-1 d-sm-none">
            <a href="https://dceg.cancer.gov/">
              <img
                src="https://analysistools.cancer.gov/common/images/DCEG-logo.svg"
                height="80"
                alt="National Cancer Institute Logo"
              />
            </a>
          </div>
        </div>
      </div>
    </header>
  );
}
