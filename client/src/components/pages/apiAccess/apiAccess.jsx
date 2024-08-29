import React from 'react';
import SwaggerUI from 'swagger-ui-react';
import './styles.scss';

export default function APIAccess() {
  return (
    <div className="mx-3">
      <div className="bg-white border p-3 mx-3">
        <div className="mb-4">
          <h1 className="h2 text-api">API Access</h1>
          <hr />
          <SwaggerUI
            url={process.env.PUBLIC_URL + '/api'}
            tryItOutEnabled={true}
          />
        </div>
      </div>
    </div>
  );
}
