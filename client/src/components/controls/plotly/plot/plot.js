import Plot from 'react-plotly.js';
import { downloadImage } from 'plotly.js';
import { Container, Row, Col, Button } from 'react-bootstrap';
import cloneDeep from 'lodash/cloneDeep';
import { saveAs } from 'file-saver';
import { handleSaveCSV } from '../../../controls/table/table2';

import './plot.scss';

export default function Plotly({
  data,
  layout,
  config,
  divId = 'plot',
  filename = 'plot',
  originalData,
  style = {},
  className = '',
  ...rest
}) {
  const defaultConfig = {
    displayModeBar: true,
    responsive: true,
    displaylogo: false,
    toImageButtonOptions: {
      format: 'svg',
      filename: filename || divId,
      scale: 1,
    },
  };

  return (
    <Container fluid style={{ minHeight: layout.height || 500 }}>
      <Row>
        <Col>
          <Plot
            className={`w-100 ${className}`}
            style={{ minHeight: layout.height || 500, ...style }}
            divId={divId}
            data={cloneDeep(data)}
            layout={cloneDeep(layout)}
            config={cloneDeep({ ...defaultConfig, ...config })}
            useResizeHandler
            {...rest}
          />
        </Col>
      </Row>
      <Row className="p-3 d-flex justify-content-end">
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
        {originalData && (
          <Col sm="auto">
            <Button
              variant="link"
              className="ml-auto m-0"
              onClick={() =>
                handleSaveCSV(originalData, `${filename || divId}.csv`)
              }
            >
              Download Data
            </Button>
          </Col>
        )}
        <Col sm="auto m-1">
          <Button
            variant="link"
            className="btn-2 btn-secondary border-12"
            onClick={() =>
              saveAs(
                new Blob([JSON.stringify({ traces: data, layout })], {
                  type: 'application/json',
                }),
                `${filename}.json`
              )
            }
          >
            Download Plotly Data &gt;
          </Button>
        </Col>
      </Row>
    </Container>
  );
}
