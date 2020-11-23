import React, { useState } from 'react';
import { Button } from 'react-bootstrap';
import { dispatchError } from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { TransformWrapper, TransformComponent } from 'react-zoom-pan-pinch';
import './plot.scss';

export default function ({ plotName, plotURL, txtPath, alt, maxHeight }) {
  const [loading, setLoading] = useState(false);

  //   download text results files
  async function downloadData(txtPath) {
    setLoading(true);
    try {
      const response = await fetch(`api/results/${txtPath}`);

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
  const zoomProps = {
    wheel: { wheelEnabled: false, step: 3 },
    zoomIn: { step: 3 },
    zoomOut: { step: 3 },
  };

  return (
    <div className="my-3">
      <LoadingOverlay active={loading} />
      <div className="d-flex align-items-end">
        <span>{plotName.slice(0, -4)}</span>
        <a className="ml-auto" href={plotURL} download={plotName}>
          Download Plot
        </a>
        {txtPath && (
          <Button
            className="p-0 border-0 ml-3"
            variant="link"
            onClick={() => downloadData(txtPath)}
          >
            Download Data
          </Button>
        )}
      </div>
      <div className="p-2 border rounded" title="Ctrl + Mouse Wheel to zoom">
        <TransformWrapper className="w-100" {...zoomProps}>
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
                  alt={alt || plotName}
                />
              </TransformComponent>
            </React.Fragment>
          )}
        </TransformWrapper>
      </div>
    </div>
  );
}
