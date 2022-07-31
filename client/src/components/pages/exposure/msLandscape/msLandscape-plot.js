import { useState, useEffect } from 'react';
import Plot from 'react-plotly.js';
import { Button, Container, Row, Col } from 'react-bootstrap';
import { cloneDeep } from 'lodash';
import { useSelector } from 'react-redux';
import { useMsLandscapePlotQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

import './plot.scss';
export default function MsLandscapePlot() {
  const publicForm = useSelector((state) => state.exposure.publicForm);
  const [params, setParams] = useState('');
  const { data, error, isFetching } = useMsLandscapePlotQuery(params, {
    skip: !params,
  });

  useEffect(() => {
    const { study, strategy, signatureSetName } = publicForm;
    if (study) {
      setParams({
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
      });
    }
  }, [publicForm]);

  return (
    <>
      <LoadingOverlay active={isFetching} />
      {data && (
        <Container fluid style={{ minHeight: '500px' }} className="mb-3">
          <Row>
            <Col>
              <Plot
                className="w-100"
                data={cloneDeep(data.traces)}
                layout={cloneDeep(data.layout)}
                config={cloneDeep(data.config)}
                useResizeHandler
              />
            </Col>
          </Row>
        </Container>
      )}
      {error && <p>An error has occured</p>}
    </>
  );
}
