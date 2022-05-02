import React from 'react';

export function Footer() {

  return (
    <footer class="py-2" style={{backgroundColor: '#003C6C'}}>
      <div class="mx-5 pl-2 pr-4">
        <div class="mt-4 text-light text-center">
          <div class="footer-nav text-left">
            <div class="mb-4 row justify-content-between">
              <div class="col-md-auto footer-nav-col">
                <div class="d-none d-lg-block footer-header text-left mb-2">
                  <h1>
                    Division of Cancer Epidemiology and Genetics
                    <span>at the National Cancer Institute</span>
                  </h1>
                </div>
                <div class="d-lg-none footer-header-mobile text-left mb-2">
                  <h1>
                    Division of Cancer Epidemiology and Genetics
                    <span>at the National Cancer Institute</span>
                  </h1>
                </div>
              </div>
            </div>
            <div class="row justify-content-between">
              <div class="col-md-auto footer-nav-col">
                <h2 class="mb-2">CONTACT INFORMATION</h2>
                <div class="my-2">
                  <a
                    class="footer-link text-light"
                    href="mailto:NCImSigPortalWebAdmin@mail.nih.gov"
                    target="_blank"
                    >Contact Us</a
                  >
                </div>
              </div>
              <div class="col-md-auto footer-nav-col">
                <h2 class="mb-2">MORE INFORMATION</h2>
                <div class="my-2">
                  <a
                    class="footer-link text-light"
                    href="https://dceg.cancer.gov/"
                    target="_blank"
                    >DCEG</a
                  >
                </div>
                <div class="my-2">
                  <a
                    class="footer-link text-light"
                    href="https://prevention.cancer.gov/"
                    target="_blank"
                    >DCP</a
                  >
                </div>
                {/* <!-- <div class="my-2">
                <a
                  class="footer-link text-light"
                  href="https://dceg.cancer.gov/research/who-we-study/cohorts/prostate-lung-colon-ovary-prospective-study"
                  target="_blank"
                  >PLCO Website</a
                >
              </div> --> */}
              </div>
              <div class="col-md-auto footer-nav-col">
                <h2 class="mb-2">POLICIES</h2>
                <div class="my-2">
                  <a
                    class="footer-link text-light"
                    href="https://www.cancer.gov/policies/accessibility"
                    target="_blank"
                    >Accessibility</a
                  >
                </div>
                <div class="my-2">
                  <a
                    class="footer-link text-light"
                    href="https://www.cancer.gov/policies/disclaimer"
                    target="_blank"
                    >Disclaimer</a
                  >
                </div>
                <div class="my-2">
                  <a
                    class="footer-link text-light"
                    href="https://www.cancer.gov/policies/foia"
                    target="_blank"
                    >FOIA</a
                  >
                </div>
                <div class="my-2">
                  <a
                    class="footer-link text-light"
                    href="https://www.cancer.gov/policies/comments"
                    target="_blank"
                    >Comment Policy</a
                  >
                </div>
              </div>
            </div>
          </div>

          <div class="footer-agencies my-4">
            <div class="row justify-content-md-center">
              <div class="col-sm-auto">
                <a
                  class="footer-link text-light"
                  href="http://www.hhs.gov/"
                  target="_blank"
                  >U.S. Department of Health and Human Services</a
                >
              </div>
              <div class="col-sm-auto">
                <a
                  class="footer-link text-light"
                  href="http://www.nih.gov"
                  target="_blank"
                  >National Institutes of Health</a
                >
              </div>
              <div class="col-sm-auto">
                <a
                  class="footer-link text-light"
                  href="https://www.cancer.gov/"
                  target="_blank"
                  >National Cancer Institute</a
                >
              </div>
              <div class="col-sm-auto">
                <a
                  class="footer-link text-light"
                  href="http://usa.gov"
                  target="_blank"
                  >USA.gov</a
                >
              </div>
            </div>
          </div>

          <div class="footer-tagline">
            <h3>
              NIH ... Turning Discovery Into Health
              <sup>®</sup>
            </h3>
          </div>
        </div>
      </div>
    </footer>
  );
}
