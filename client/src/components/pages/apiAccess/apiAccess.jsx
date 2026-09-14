import SwaggerUI from 'swagger-ui-react';
import SwaggerLabelInjector from './swagger-ui/swagger-injecttion.js';
import SwaggerColorCustomizer from './swagger-ui/swagger-color-customizer.js';
import SwaggerScrollablePreEnhancer from './swagger-ui//swagger-scrollable.js';

import './styles.scss';

export default function APIAccess() {
  return (
    <div className="mx-3">
      <div className="bg-white border p-3 mx-3">
        <div className="mb-4">
          <h1 className="h2 text-api">API Access</h1>
          <p>
            This page provides an interactive interface to the mSigPortal REST
            API. Browse the available endpoints below, expand one to see its
            parameters, then select "Try it out" to send a live request and view
            the response.
          </p>
          <p className="text-muted small mb-0">
            Requests are rate limited to 3,000 requests per 15 minutes per IP
            address. Exceeding this limit returns an HTTP 429 response;
            remaining quota is reported in the standard <code>RateLimit</code>{' '}
            response headers.
          </p>
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
