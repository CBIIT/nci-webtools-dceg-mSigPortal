import React from 'react';
import SwaggerUI from 'swagger-ui-react';

export default function APIAccess() {
  return (
    <div className="mx-3">
      <div className="bg-white border p-3 mx-3">
        <div className="mb-4">
          <h1 class="h2">API Access</h1>
          <hr />
          <SwaggerUI url={process.env.PUBLIC_URL + '/api'} />
        </div>
      </div>
    </div>
  );
}
