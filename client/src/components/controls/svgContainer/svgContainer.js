import React, { useState } from 'react';
import { Button } from 'react-bootstrap';
import { LoadingOverlay } from '../loading-overlay/loading-overlay';
import { TransformWrapper, TransformComponent } from 'react-zoom-pan-pinch';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faPlus, faMinus, faHome } from '@fortawesome/free-solid-svg-icons';
import { useDispatch } from 'react-redux';
import { actions } from '../../../services/store/modal';
import './styles.scss';

export default function SvgContainer({
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
    pinch: { step: 0.1 },
    wheel: { step: 0.01, activationKeys: ['z'] },
  };

  // fetch image to refresh cached image in chromium browsers
  // need to do this becuase plotpaths are always the same for recalculations
  try {
    fetch(plotPath, { cache: 'reload', mode: 'no-cors' });
  } catch (_) {}

  return (
    <div
      className={`${className} mx-auto`}
      title="Z + Mouse Wheel to zoom"
      style={{ width: 'auto', height: '100%', maxWidth: '1500px' }}
    >
      <LoadingOverlay active={loading} />
      <TransformWrapper {...zoomProps}>
        {({ zoomIn, zoomOut, resetTransform }) => (
          <React.Fragment>
            {plotPath && (
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
                    onClick={() => zoomIn(0.2)}
                    aria-label="zoom in"
                  >
                    <FontAwesomeIcon
                      icon={faPlus}
                      style={{ color: '#fafafa' }}
                    />
                  </Button>
                  <Button
                    size="sm"
                    className="ml-1"
                    variant="secondary"
                    onClick={() => zoomOut(0.2)}
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
                    onClick={() => resetTransform()}
                    aria-label="reset zoom"
                  >
                    <FontAwesomeIcon
                      icon={faHome}
                      style={{ color: '#fafafa' }}
                    />
                  </Button>
                </div>
              </div>
            )}
            <TransformComponent>
              <img
                className="w-100"
                src={
                  plotPath + (cacheBreaker ? `#${new Date().getTime()}` : '')
                }
                style={{ maxHeight: height || '600px' }}
                alt={alt || 'Plot Unavailable'}
              />
            </TransformComponent>
          </React.Fragment>
        )}
      </TransformWrapper>
    </div>
  );
}
