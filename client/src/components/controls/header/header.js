import React from 'react';
import { useLocation } from 'react-router-dom';

export function Header() {
  const location = useLocation();
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
        <div
          className={
            location.pathname === '/'
              ? 'container d-none d-lg-block'
              : 'd-none d-lg-block'
          }
        >
          <div className="">
            <div className="">
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

              <div className="p-1 d-block d-sm-none">
                <a href="https://dceg.cancer.gov/">
                  <img
                    src="https://analysistools.cancer.gov/common/images/DCEG-logo.svg"
                    height="80"
                    width="100%"
                    alt="National Cancer Institute Logo"
                  />
                </a>
              </div>
            </div>
          </div>
        </div>

        <div className=" d-block d-lg-none">
          <div className="">
            <div className="">
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

              <div className="p-1 d-block d-sm-none">
                <a href="https://dceg.cancer.gov/">
                  <img
                    src="https://analysistools.cancer.gov/common/images/DCEG-logo.svg"
                    height="80"
                    width="100%"
                    alt="National Cancer Institute Logo"
                  />
                </a>
              </div>
            </div>
          </div>
        </div>
      </div>
    </header>
  );
}
