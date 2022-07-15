import { useState, useEffect } from 'react';
import Plot from 'react-plotly.js';
import { Button, Container, Row, Col } from 'react-bootstrap';
import { cloneDeep } from 'lodash';
import { useSelector } from 'react-redux';
import { useTmbPlotQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

import './plot.scss';
export default function MutProfilePlot() {
  const publicForm = useSelector((state) => state.exposure.publicForm);
  const [params, setParams] = useState('');
  const { data, error, isFetching } = useTmbPlotQuery(params, {
    skip: !params,
  });

  useEffect(() => {
    const { study, strategy, signatureSetName, signatureName } = publicForm;
    if (study) {
      setParams({
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        signatureName: signatureName.value,
      });
    }
  }, [publicForm]);
  //console.log(data);
  //console.log(publicForm);
  if (data) {
    //console.log(data.traces.length);
  }

  return (
    <>
      <LoadingOverlay active={isFetching} />
      {data && (
        <Container fluid style={{ minHeight: '500px' }} className="mb-3">
          <Row className="justify-content-center text-center">
            <Col
              {...(data.traces.length > 1
                ? { className: '' }
                : { className: 'col-md-3 content col-md-offset-3 ' })}
            >
              <Plot
                // {...(data.traces.length > 1
                //   ? { className: 'w-100' }
                //   : { className: 'w-30' })}
                // style={{ height: '500px' }}
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
