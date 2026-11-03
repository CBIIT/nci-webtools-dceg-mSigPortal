import { useState } from 'react';
import { Accordion, Card } from 'react-bootstrap';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faChevronDown, faChevronUp } from '@fortawesome/free-solid-svg-icons';
import SwaggerUI from 'swagger-ui-react';
import SwaggerLabelInjector from './swagger-ui/swagger-injecttion.js';
import SwaggerColorCustomizer from './swagger-ui/swagger-color-customizer.js';
import SwaggerScrollablePreEnhancer from './swagger-ui//swagger-scrollable.js';

import './styles.scss';

export default function APIAccess() {
  const [open, setOpen] = useState(false);

  const rExample = `library(httr)
library(jsonlite)

# Retrieve mutational spectrum options for a study
res <- GET(
  "https://analysistools.cancer.gov/mutational-signatures/api/mutational_spectrum_options",
  query = list(study = "PCAWG", cancer = "Lung-AdenoCA", profile = "SBS")
)

data <- fromJSON(content(res, "text", encoding = "UTF-8"))
print(data)`;

  const pythonExample = `import requests

# Retrieve mutational spectrum options for a study
res = requests.get(
    "https://analysistools.cancer.gov/mutational-signatures/api/mutational_spectrum_options",
    params={"study": "PCAWG", "cancer": "Lung-AdenoCA", "profile": "SBS"},
)

data = res.json()
print(data)`;

  return (
    <div className="mx-3">
      <div className="bg-white border p-3 mx-3">
        <div className="mb-4">
          <h1 className="h2 text-api">API Access</h1>
          <div id="description">
            <p>
              This page provides an interactive interface to the mSigPortal REST
              API. Browse the available endpoints below, expand one to see its
              parameters, then select "Try it out" to send a live request and
              view the response.
            </p>
            <p className="mb-0">
              Requests are rate limited to 3,000 requests per 15 minutes per IP
              address. Exceeding this limit returns an HTTP 429 response;
              remaining quota is reported in the standard <code>RateLimit</code>{' '}
              response headers.
            </p>
          </div>
          <Accordion className="mt-3">
            <Card>
              <Accordion.Toggle
                as={Card.Header}
                eventKey="0"
                className="font-weight-bold d-flex justify-content-between"
                style={{ cursor: 'pointer' }}
                onClick={() => setOpen(!open)}
              >
                Usage Example
                <FontAwesomeIcon icon={open ? faChevronUp : faChevronDown} />
              </Accordion.Toggle>
              <Accordion.Collapse eventKey="0">
                <Card.Body>
                  <p className="small text-muted">
                    Basic examples calling the{' '}
                    <code>/api/mutational_spectrum_options</code> endpoint.
                  </p>
                  <h6 className="font-weight-bold">R</h6>
                  <pre
                    className="bg-light border rounded p-3 mb-3 small"
                    style={{ whiteSpace: 'pre', overflowX: 'auto' }}
                  >
                    <code>{rExample}</code>
                  </pre>
                  <h6 className="font-weight-bold">Python</h6>
                  <pre
                    className="bg-light border rounded p-3 mb-0 small"
                    style={{ whiteSpace: 'pre', overflowX: 'auto' }}
                  >
                    <code>{pythonExample}</code>
                  </pre>
                </Card.Body>
              </Accordion.Collapse>
            </Card>
          </Accordion>
          <hr />
          <SwaggerUI
            url={import.meta.env.BASE_URL + '/api'}
            tryItOutEnabled={true}
          />
          {/* Add SwaggerLabelInjector to observe and inject the label */}
          <SwaggerLabelInjector />
          {/* Custom integer color styling */}
          <SwaggerColorCustomizer />
          <SwaggerScrollablePreEnhancer />
        </div>
      </div>
    </div>
  );
}
