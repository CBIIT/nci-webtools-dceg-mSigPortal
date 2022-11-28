import Plot from 'react-plotly.js';
import { downloadImage } from 'plotly.js';
import { Container, Row, Col, Button } from 'react-bootstrap';
import cloneDeep from 'lodash/cloneDeep';
import { saveAs } from 'file-saver';

import './plot.scss';

export default function Plotly({
  data,
  layout,
  config,
  divId,
  filename,
  originalData,
  ...rest
}) {
  const defaultConfig = {
    displayModeBar: true,
    responsive: true,
    displaylogo: false,
    toImageButtonOptions: {
      format: 'svg',
      filename: filename || divId || 'plot',
      scale: 1,
    },
  };

  return (
    <Container fluid>
      <Row>
        {/* <Col sm="auto">
          <Button
            variant="link"
            onClick={() =>
              downloadImage(divId, {
                format: 'png',
                filename: filename,
              })
            }
          >
            Download PNG
          </Button>
        </Col>
        <Col sm="auto">
          <Button
            variant="link"
            onClick={() =>
              downloadImage(divId, {
                format: 'svg',
                filename: filename,
              })
            }
          >
            Download SVG
          </Button>
        </Col> */}
        <Col className="d-flex justify-content-end">
          <Button
            variant="link"
            onClick={() =>
              saveAs(
                new Blob([JSON.stringify({ traces: data, layout })], {
                  type: 'application/json',
                }),
                `${filename}.json`
              )
            }
          >
            Download Plotly Data
          </Button>
        </Col>
      </Row>
      <Row>
        <Col>
          <Plot
            className="w-100"
            divId={divId}
            data={cloneDeep(data)}
            layout={cloneDeep(layout)}
            config={cloneDeep({ ...defaultConfig, ...config })}
            useResizeHandler
            {...rest}
          />
        </Col>
      </Row>
    </Container>
  );
}
