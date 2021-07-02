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
  maxHeight,
  className,
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
    wheel: { wheelEnabled: false, step: 3 },
    zoomIn: { step: 3 },
    zoomOut: { step: 3 },
  };

  return (
    <div>
      <LoadingOverlay active={loading} />
      <div className={className} title="Ctrl + Mouse Wheel to zoom">
        <TransformWrapper className="w-100" {...zoomProps}>
          {({ zoomIn, zoomOut, resetTransform }) => (
            <React.Fragment>
              <Row className="tools">
                <Col />
                <Col className="d-flex text-nowrap">
                  {title && <strong className="mx-auto">{title}</strong>}
                </Col>
                <Col className="d-flex text-nowrap">
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
                        onClick={() => download(`api/results/${txtPath}`)}
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
                    <FontAwesomeIcon
                      icon={faPlus}
                      style={{ color: '#fafafa' }}
                    />
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
                    <FontAwesomeIcon
                      icon={faHome}
                      style={{ color: '#fafafa' }}
                    />
                  </Button>
                </Col>
              </Row>

              <TransformComponent>
                <img
                  className="w-100"
                  src={`${plotPath}#${Date.now()}`}
                  style={{ maxHeight: maxHeight || '500px' }}
                  alt={alt || 'Plot Unavailable'}
                />
              </TransformComponent>
            </React.Fragment>
          )}
        </TransformWrapper>
      </div>
    </div>
  );
}
