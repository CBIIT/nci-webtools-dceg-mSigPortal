import React from 'react';
import SwaggerUI from 'swagger-ui-react';
import SwaggerLabelInjector from "./swagger-ui/swagger-injecttion.js";
import SwaggerColorCustomizer from "./swagger-ui/swagger-color-customizer.js";
import SwaggerScrollablePreEnhancer from "./swagger-ui//swagger-scrollable.js";

import './styles.scss';

export default function APIAccess() {
  return (
    <div className="mx-3">
      <div className="bg-white border p-3 mx-3">
        <div className="mb-4">
          <h1 className="h2 text-api">API Access</h1>
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
