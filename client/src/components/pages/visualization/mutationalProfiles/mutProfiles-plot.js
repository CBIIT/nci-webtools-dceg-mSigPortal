import { useEffect, useState } from 'react';
import cloneDeep from 'lodash/cloneDeep';
import { Button, Container, Row, Col } from 'react-bootstrap';
import Plot from 'react-plotly.js';
import { downloadImage } from 'plotly.js';
import { saveAs } from 'file-saver';
import { useSelector } from 'react-redux';
import { useMutationalProfilesQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

import './plot.scss';

export default function MutProfilePlot() {
  const store = useSelector((state) => state.visualization);

  const { sample, profile, matrix } = store.mutationalProfiles;
  const { study, cancer, strategy } = store.publicForm;
  const { source, projectID } = store.main;

  const [params, setParams] = useState(null);

  const { data, error, isFetching } = useMutationalProfilesQuery(params, {
    skip: !params,
  });
  console.log(store.main);

  // get data on form change
  useEffect(() => {
    if (sample && sample.value && profile.value && matrix.value) {
      const params = {
        study: study.value,
        cancer: cancer.value,
        strategy: strategy.value,
        ...(source == 'user' && { userId: projectID }),
        sample: sample.value,
        profile: profile.value,
        matrix: matrix.value,
      };
      setParams(params);
    }
  }, [sample, profile, matrix]);

  const divId = 'mutationalProfilePlot';
  const config = {
    displayModeBar: true,
    responsive: true,
    displaylogo: false,
    toImageButtonOptions: {
      format: 'svg',
      filename: sample?.value || 'Mutational Profile',
      scale: 1,
    },
  };

  return (
    <>
      <LoadingOverlay active={isFetching} />
      {error ? (
        <div className="text-center">
          <div>An error as occured</div>
          <div>{error.message}</div>
        </div>
      ) : (
        data && (
          <Container fluid style={{ minHeight: '500px' }} className="mb-3">
            <Row>
              <Col>
                <Plot
                  className="w-100"
                  divId={divId}
                  data={cloneDeep(data.traces)}
                  layout={cloneDeep(data.layout)}
                  config={cloneDeep(config)}
                  useResizeHandler
                />
              </Col>
            </Row>
            {/* <Row className="justify-content-center"> */}
            {/* <Col sm="auto">
                <Button
                  onClick={() =>
                    downloadImage(divId, {
                      format: "png",
                      filename: sample.value,
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
                      format: "svg",
                      filename: sample.value,
                    })
                  }
                >
                  Download SVG
                </Button>
              </Col> */}
            <Row>
              <Col sm="auto">
                <Button
                  onClick={() =>
                    saveAs(
                      new Blob([JSON.stringify(data)], {
                        type: 'application/json',
                      }),
                      `${sample.value}.json`
                    )
                  }
                >
                  Download JSON
                </Button>
              </Col>
            </Row>
          </Container>
        )
      )}
    </>
  );
}
