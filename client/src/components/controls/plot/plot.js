import React, { useState } from 'react';
import { Button } from 'react-bootstrap';
import { dispatchError } from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { TransformWrapper, TransformComponent } from 'react-zoom-pan-pinch';

export default function ({ plotName, plotURL, txtPath, maxHeight }) {
  const rootURL = window.location.pathname;

  const [loading, setLoading] = useState(false);

  //   download text results files
  async function downloadResults(txtPath) {
    setLoading(true);
    try {
      const response = await fetch(`${rootURL}downloadPlotData`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ path: txtPath }),
      });
      if (!response.ok) {
        const { msg } = await response.json();
        dispatchError(msg);
      } else {
        const file = await response.blob();
        const objectURL = URL.createObjectURL(file);
        const tempLink = document.createElement('a');

        tempLink.href = objectURL;
        tempLink.download = txtPath.split('/').slice(-1)[0];
        document.body.appendChild(tempLink);
        tempLink.click();
        document.body.removeChild(tempLink);
        URL.revokeObjectURL(objectURL);
      }
    } catch (err) {
      dispatchError(err);
    }
    setLoading(false);
  }

  return (
    <div>
      <LoadingOverlay active={loading} />
      <div className="d-flex">
        <a className="px-2 py-1" href={plotURL} download={plotName}>
          Download Plot
        </a>
        {txtPath && (
          <span className="ml-auto">
            <Button
              className="px-2 py-1"
              variant="link"
              onClick={() => downloadResults(txtPath)}
            >
              Download Results
            </Button>
          </span>
        )}
      </div>
      <div className="p-2 border rounded">
        <TransformWrapper defaultScale={1}>
          {({ zoomIn, zoomOut, resetTransform }) => (
            <React.Fragment>
              <div className="tools">
                <Button variant="secondary" onClick={zoomIn}>
                  +
                </Button>{' '}
                <Button variant="secondary" onClick={zoomOut}>
                  -
                </Button>{' '}
                <Button variant="secondary" onClick={resetTransform}>
                  Reset
                </Button>
              </div>
              <TransformComponent>
                <img
                  className="w-100"
                  src={plotURL}
                  style={{ maxHeight: maxHeight || '500px' }}
                />
              </TransformComponent>
            </React.Fragment>
          )}
        </TransformWrapper>
      </div>
    </div>
  );
}
