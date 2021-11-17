import React from 'react';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

export default function Instructions({ loading }) {
  return (
    <div className="bg-white border rounded py-3 px-4">
      <LoadingOverlay active={loading} />
      <h4>Instructions</h4>
      <p>
        Choose a Data Source and its associated options to submit a query using
        the panel on the left
      </p>
      <hr />
      <h4>Data Source</h4>
      <p>Public: Perform analysis using data available on the website</p>
      {/* <p>User: Upload your own data</p> */}
    </div>
  );
}
