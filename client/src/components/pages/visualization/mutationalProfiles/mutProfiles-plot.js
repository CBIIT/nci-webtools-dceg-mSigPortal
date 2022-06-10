import cloneDeep from 'lodash/cloneDeep';
import { useRecoilValue } from 'recoil';
import { Button, Row, Col } from 'react-bootstrap';
import { getPlot, formState } from './mutProfiles.state';
import Plot from 'react-plotly.js';
import { downloadImage } from 'plotly.js';
import { saveAs } from 'file-saver';

export default function MutProfilePlot() {
  const { data, layout, config } = useRecoilValue(getPlot);
  const { option } = useRecoilValue(formState);

  const divId = 'mutationalProfilePlot';
  return (
    <div>
      {data.length ? (
        <div>
          <Plot
            className="w-100"
            divId={divId}
            style={{ height: '500px' }}
            data={cloneDeep(data)}
            layout={cloneDeep(layout)}
            config={cloneDeep(config)}
            useResizeHandler
          />
          <Row className="justify-content-center">
            <Col sm="auto">
              <Button
                onClick={() =>
                  downloadImage(divId, {
                    format: 'png',
                    filename: option.value.signature,
                  })
                }
              >
                Download PNG
              </Button>
            </Col>
            <Col sm="auto">
              <Button
                onClick={() =>
                  downloadImage(divId, {
                    format: 'svg',
                    filename: option.value.signature,
                  })
                }
              >
                Download SVG
              </Button>
            </Col>
            <Col sm="auto">
              <Button
                onClick={() =>
                  saveAs(
                    new Blob([JSON.stringify(data)], {
                      type: 'application/json',
                    }),
                    `${option.value.signature}.json`
                  )
                }
              >
                Download JSON
              </Button>
            </Col>
          </Row>
        </div>
      ) : (
        <div>Select a Profile</div>
      )}
    </div>
  );
}
