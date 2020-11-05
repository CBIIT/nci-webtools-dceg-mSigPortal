import React, { useState, useEffect } from 'react';
import { Button } from 'react-bootstrap';
import { dispatchError } from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { TransformWrapper, TransformComponent } from 'react-zoom-pan-pinch';
import './plot.scss';

export default function ({ plotName, plotURL, txtPath, alt, maxHeight }) {
  const rootURL = window.location.pathname;

  const [loading, setLoading] = useState(false);

  // Key listener hook usehooks.com/useKeyPress/
  const useKeyPress = (targetKey) => {
    // State for keeping track of whether key is pressed
    const [keyPressed, setKeyPressed] = useState(false);
    // If pressed key is our target key then set to true
    const downHandler = ({ key }) => {
      if (key === targetKey) {
        console.log(key + 'DOWN');
        setKeyPressed(true);
      }
    };

    // If released key is our target key then set to false
    const upHandler = ({ key }) => {
      if (key === targetKey) {
        console.log(key + 'UP');

        setKeyPressed(false);
      }
    };

    // Add event listeners
    useEffect(() => {
      window.addEventListener('keydown', downHandler);
      window.addEventListener('keyup', upHandler);
      // Remove event listeners on cleanup
      return () => {
        window.removeEventListener('keydown', downHandler);
        window.removeEventListener('keyup', upHandler);
      };
    }, []); // Empty array ensures that effect is only run on mount and unmount

    return keyPressed;
  };

  const zoomKeyDown = useKeyPress('z');

  //   download text results files
  async function downloadData(txtPath) {
    setLoading(true);
    try {
      const response = await fetch(`${rootURL}results/${txtPath}`);

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
    wheel: { wheelEnabled: zoomKeyDown, step: 3 },
    zoomIn: { step: 3 },
    zoomOut: { step: 3 },
  };

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
              onClick={() => downloadData(txtPath)}
            >
              Download Data
            </Button>
          </span>
        )}
      </div>
      <div className="p-2 border rounded">
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
