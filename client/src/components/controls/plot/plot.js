import React, { useState } from 'react';
import { Button, Row, Col } from 'react-bootstrap';
import { dispatchError } from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { TransformWrapper, TransformComponent } from 'react-zoom-pan-pinch';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faPlus, faMinus, faHome } from '@fortawesome/free-solid-svg-icons';
import './plot.scss';

export default function Plot({
  plotName,
  plotURL,
  txtPath,
  alt,
  maxHeight,
  className,
  ...rest
}) {
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
    <div>
      <LoadingOverlay active={loading} />
      <div className={className} title="Ctrl + Mouse Wheel to zoom">
        <TransformWrapper className="w-100" {...zoomProps}>
          {({ zoomIn, zoomOut, resetTransform }) => (
            <React.Fragment>
              <Row className="tools">
                <Col />
                <Col className="d-flex">
                  {plotName && <strong className="mx-auto">{plotName}</strong>}
                </Col>
                <Col className="d-flex">
                  <div className="d-flex align-items-end ml-auto mb-auto">
                    <a
                      className="ml-auto"
                      href={plotURL}
                      download={plotName.trim()}
                    >
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
                  <Button
                    size="sm"
                    className="ml-3"
                    variant="secondary"
                    onClick={zoomIn}
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
