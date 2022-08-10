import { useEffect, useState } from 'react';
import cloneDeep from 'lodash/cloneDeep';
import { NavHashLink } from 'react-router-hash-link';
import { Button, Container, Row, Col } from 'react-bootstrap';
import Plot from 'react-plotly.js';
import { useSelector } from 'react-redux';
import {
  usePcBetweenSamplesQuery,
  usePcToReferenceSignaturesQuery,
} from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import './plot.scss';

export default function ProfileComparisonPlot() {
  const publicForm = useSelector((state) => state.visualization.publicForm);
  const [params, setParams] = useState('');
  const store = useSelector((state) => state.visualization);
  const { study, cancer, strategy } = store.publicForm;
  const { profile, sample1, sample2 } = store.profileComparison;
  const { source } = store.main;

  const { data, error, isFetching } = usePcBetweenSamplesQuery(params, {
    skip: !params,
  });

  useEffect(() => {
    if (sample1 && sample1.value && sample2 && sample2.value && profile.value) {
      const params = {
        study: study.value,
        cancer: cancer.value,
        strategy: strategy.value,
        sample1: sample1.value,
        sample2: sample2.value,
        profile: profile.value,
      };
      setParams(params);
    }
  }, [sample1, sample2, profile]);

  const divId = 'pcBetweenSamples';
  const config = {
    displayModeBar: true,
    responsive: true,
    displaylogo: false,
    toImageButtonOptions: {
      format: 'svg',
      filename: profile?.value || 'Mutational Profile',
      scale: 1,
    },
  };

  return (
    <>
      <LoadingOverlay active={isFetching} />
      {error && (
        <p className="p-3">An error has occured. Please verify your input.</p>
      )}
      <div id="barchart">
        {data && (
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
            <Row>
              <div className="p-3">
                <p>
                  The plot above shows the mutational profiles of two selected
                  samples, as well as the difference between them. The text at
                  the top of the plot indicates the profile similarity
                  calculated using Residual Sum of Squares (RSS) and cosine
                  similarity methods.
                </p>
                <p>
                  RSS measures the discrepancy between two mutational profiles.
                  Cosine similarity measures how similar two mutational profiles
                  are. For example, two identical mutational profiles will have
                  RSS = 0 and Cosine similarity = 1. For additional information
                  about RSS and cosine similarity, click{' '}
                  <NavHashLink to="/faq#cosine-similarity">here</NavHashLink>.
                </p>
              </div>
            </Row>
          </Container>
        )}
        {/* {patternData?.output.plotPath && !patternData?.output.barPath && (
          <div className="p-3">
            <p>Frequency of Mutational Pattern</p>
            <p>
              No mutational pattern with proportion of mutations large than{' '}
              {proportion}
            </p>
          </div>
        )} */}
      </div>{' '}
    </>
  );
}
