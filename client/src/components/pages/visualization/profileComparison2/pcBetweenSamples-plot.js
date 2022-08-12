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
  const [params1, setParams1] = useState('');
  const [params2, setParams2] = useState('');
  const store = useSelector((state) => state.visualization);
  const { study, cancer, strategy } = store.publicForm;
  const { profile, sample1, sample2 } = store.profileComparison;
  const { source } = store.main;

  const {
    data: sample1data,
    error: sample1error,
    isFetching: sample1fetching,
  } = usePcBetweenSamplesQuery(params1, {
    skip: !params1,
  });

  const {
    data: sample2data,
    error: sample2error,
    isFetching: sample2fetching,
  } = usePcBetweenSamplesQuery(params2, {
    skip: !params2,
  });

  useEffect(() => {
    if (sample1 && sample1.value && profile.value) {
      const params1 = {
        study: study.value,
        cancer: cancer.value,
        strategy: strategy.value,
        sample: sample1.value,
        profile: profile.value,
      };
      setParams1(params1);
    }
  }, [sample1, profile]);

  useEffect(() => {
    if (sample2 && sample2.value && profile.value) {
      const params2 = {
        study: study.value,
        cancer: cancer.value,
        strategy: strategy.value,
        sample: sample2.value,
        profile: profile.value,
      };
      setParams2(params2);
    }
  }, [sample2, profile]);

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
      <LoadingOverlay active={sample1fetching} />
      {sample1error && (
        <p className="p-3">An error has occured. Please verify your input.</p>
      )}
      <div id="barchart">
        {sample1data &&
          sample2data(
            <Container fluid style={{ minHeight: '500px' }} className="mb-3">
              <Row>
                <Col>
                  <Plot
                    className="w-100"
                    divId={divId}
                    data={cloneDeep([sample1data, sample2data].traces)}
                    layout={cloneDeep([sample1data, sample2data].layout)}
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
                    RSS measures the discrepancy between two mutational
                    profiles. Cosine similarity measures how similar two
                    mutational profiles are. For example, two identical
                    mutational profiles will have RSS = 0 and Cosine similarity
                    = 1. For additional information about RSS and cosine
                    similarity, click{' '}
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
