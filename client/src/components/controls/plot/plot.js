import React, { useState } from 'react';
import { Button, Row, Col } from 'react-bootstrap';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { TransformWrapper, TransformComponent } from 'react-zoom-pan-pinch';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faPlus, faMinus, faHome } from '@fortawesome/free-solid-svg-icons';
import { useDispatch } from 'react-redux';
import { actions } from '../../../services/store/modal';
import './plot.scss';

export default function Plot({
  title,
  downloadName,
  plotPath,
  txtPath,
  alt,
  height,
  className,
  cacheBreaker = true,
  ...rest
}) {
  const dispatch = useDispatch();
  const mergeError = (state) => dispatch(actions.mergeModal({ error: state }));

  const [loading, setLoading] = useState(false);

  async function download(path) {
    setLoading(true);
    try {
      const response = await fetch(path);

      if (!response.ok) {
        const { msg } = await response.json();
        mergeError({ visible: true, message: msg });
      } else {
        const file = await response.blob();
        const objectURL = URL.createObjectURL(file);
        const tempLink = document.createElement('a');

        tempLink.href = objectURL;
        tempLink.download = path.split('/').slice(-1)[0];
        document.body.appendChild(tempLink);
        tempLink.click();
        document.body.removeChild(tempLink);
        URL.revokeObjectURL(objectURL);
      }
    } catch (err) {
      mergeError({ visible: true, message: err.message });
    }
    setLoading(false);
  }

  const zoomProps = {
    wheel: { wheelEnabled: false, step: 5 },
    zoomIn: { step: 5 },
    zoomOut: { step: 5 },
  };

  return (
    <div
      className={`${className}`}
      title="Ctrl + Mouse Wheel to zoom"
      style={{ width: 'auto', height: '100%' }}
    >
      <LoadingOverlay active={loading} />
      <TransformWrapper {...zoomProps}>
        {({ zoomIn, zoomOut, resetTransform }) => (
          <React.Fragment>
            <div className="tools mb-3">
              <div className="d-flex justify-content-center">
                {title && <strong className="mb-3">{title}</strong>}
              </div>
              <div className="d-flex">
                <div className="d-flex align-items-end ml-auto mb-auto">
                  <Button
                    className="p-0 border-0 ml-3"
                    variant="link"
                    onClick={() => download(plotPath)}
                  >
                    Download Plot
                  </Button>
                  {txtPath && (
                    <Button
                      className="p-0 border-0 ml-3"
                      variant="link"
                      onClick={() => download(txtPath)}
                    >
                      Download Data
                    </Button>
                  )}
                </div>
                <Button
                  size="sm"
                  className="ml-3"
                  variant="secondary"
                  onClick={zoomIn}
                  aria-label="zoom in"
                >
                  <FontAwesomeIcon icon={faPlus} style={{ color: '#fafafa' }} />
                </Button>
                <Button
                  size="sm"
                  className="ml-1"
                  variant="secondary"
                  onClick={zoomOut}
                  aria-label="zoom out"
                >
                  <FontAwesomeIcon
                    icon={faMinus}
                    style={{ color: '#fafafa' }}
                  />
                </Button>
                <Button
                  size="sm"
                  className="ml-1"
                  variant="secondary"
                  onClick={resetTransform}
                  aria-label="reset zoom"
                >
                  <FontAwesomeIcon icon={faHome} style={{ color: '#fafafa' }} />
                </Button>
              </div>
            </div>

            <TransformComponent>
              <img
                className="w-100"
                src={
                  plotPath + `${cacheBreaker ? `#${new Date().getTime()}` : ''}`
                }
                style={{ maxHeight: height || '500px' }}
                alt={alt || 'Plot Unavailable'}
              />
            </TransformComponent>
          </React.Fragment>
        )}
      </TransformWrapper>
    </div>
  );
}
