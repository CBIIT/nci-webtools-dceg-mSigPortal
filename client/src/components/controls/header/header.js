import React from 'react';

export function Header() {

  return (
    <header class="bg-light">
      <div>
        <a
          href="/#/#root"
          class="
            sr-only sr-only-focusable
            d-block
            text-white
            bg-primary-dark
            text-center
          "
          >Skip to Main Content</a
        >
        <div class="ml-4 mr-5">
          <div class="p-2 d-none d-sm-block">
            <a href="https://dceg.cancer.gov/" target="_blank">
              <img
                src="https://analysistools.cancer.gov/common/images/DCEG-logo.svg"
                height="100"
                alt="National Cancer Institute Logo"
              />
            </a>
            {/* <a href="/#/">
              <img
                class="d-none d-md-block"
                src="assets/images/msigportal-logo.png"
                alt="mSigPortal Logo"
                height="80"
                style={{float: 'right', marginTop: '15px'}}
              />
            </a> */}
          </div>
          <div class="p-1 d-sm-none">
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
