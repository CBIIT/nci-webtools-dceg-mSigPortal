import React, { useState, useEffect } from 'react';
import { Form, Row, Col, Button, Card } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import Plot from 'react-plotly.js';
import { useSatsDataQuery, useSatsOptionsQuery } from './satsApiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Select from '../../../controls/select/selectHookForm';

export default function SATSSection() {
  const [params, setParams] = useState(null);

  // Fetch sample etiology data to get available studies and signature sets
  const { data: options, isFetching: fetchingOptions } = useSatsOptionsQuery();

  // Fetch SATS plot data
  const { 
    data: plotConfig, 
    isFetching: fetchingPlot, 
    error: plotError 
  } = useSatsDataQuery(params, { skip: !params });

  const { control, handleSubmit, watch, setValue } = useForm({
    defaultValues: {
      study: null,
      signatureSetName: null,
    }
  });

  const { study, signatureSetName } = watch();

  // Create form options from etiology data
  const studyOptions = options
    ? [...new Set(options.map(e => e.study))]
        .filter(Boolean)
        .map(study => ({ label: study, value: study }))
    : [];

  const signatureSetOptions = options && study
    ? [...new Set(options
        .filter(e => e.study === study.value)
        .map(e => e.signatureSetName))]
        .filter(Boolean)
        .map(set => ({ label: set, value: set }))
    : [];

  // Set initial values
  useEffect(() => {
    if (studyOptions.length > 0 && !study) {
      setValue('study', studyOptions[0]);
    }
  }, [studyOptions, study, setValue]);

  useEffect(() => {
    if (signatureSetOptions.length > 0 && !signatureSetName) {
      setValue('signatureSetName', signatureSetOptions[0]);
    }
  }, [signatureSetOptions, signatureSetName, setValue]);

  function onSubmit(data) {
    const queryParams = {
      study: data.study?.value,
      signatureSetName: data.signatureSetName?.value,
    };
    setParams(queryParams);
  }

  return (
    <div className="mt-4">
      <h5 className="separator">SATS - Signature Activity in Tumor Samples</h5>
      <p className="text-muted mb-3">
        Interactive visualization showing tumor mutational burden (TMB) and signature presence 
        across different cancer types.
      </p>

      <Card className="mb-3">
        <Card.Body>
          <LoadingOverlay active={fetchingOptions} />
          <Form onSubmit={handleSubmit(onSubmit)}>
            <Row>
              <Col md={4}>
                <Select
                  name="study"
                  label="Study"
                  options={studyOptions}
                  control={control}
                  disabled={fetchingOptions}
                />
              </Col>
              <Col md={4}>
                <Select
                  name="signatureSetName"
                  label="Signature Set"
                  options={signatureSetOptions}
                  control={control}
                  disabled={fetchingOptions || !study}
                />
              </Col>
              <Col md={4} className="d-flex">
                <Button
                  type="submit"
                  variant="primary"
                  className="mt-auto mb-3"
                  disabled={fetchingOptions || !study || !signatureSetName}
                >
                  Generate SATS Plot
                </Button>
              </Col>
            </Row>
          </Form>

          {params && (
            <div className="mt-3 p-2 bg-light rounded">
              <small>
                <strong>Current Selection:</strong> {params.study} - {params.signatureSetName}
              </small>
            </div>
          )}
        </Card.Body>
      </Card>

      <Card>
        <Card.Body>
          <LoadingOverlay active={fetchingPlot} />
          
          {plotError && (
            <div className="alert alert-danger">
              <h6>Error loading plot data</h6>
              <p>{plotError.message || 'An unknown error occurred'}</p>
            </div>
          )}

          {plotConfig && plotConfig.traces && plotConfig.traces.length > 0 ? (
            <Plot
              data={plotConfig.traces}
              layout={plotConfig.layout}
              config={plotConfig.config}
              style={{ width: '100%', height: '800px' }}
            />
          ) : params && !fetchingPlot && !plotError ? (
            <div className="text-center p-5 text-muted">
              <h6>No Data Available</h6>
              <p>No signature presence data found for the selected study and signature set.</p>
              <p>Try selecting a different study or signature set.</p>
            </div>
          ) : !params ? (
            <div className="text-center p-5 text-muted">
              <h6>Select Data to View Plot</h6>
              <p>Choose a study and signature set above to generate the SATS plot.</p>
            </div>
          ) : null}
        </Card.Body>
      </Card>

      {/* Info section */}
      <Card className="mt-3">
        <Card.Body>
          <Row>
            <Col md={6}>
              <h6>Plot Components:</h6>
              <ul>
                <li><strong>Top Panel:</strong> Stacked bar chart showing tumor mutational burden (TMB) per megabase for each cancer type</li>
                <li><strong>Bottom Panel:</strong> Dot plot showing the proportion of samples in each cancer type that exhibit each signature</li>
              </ul>
            </Col>
            <Col md={6}>
              <h6>Interpretation:</h6>
              <ul>
                <li><strong>Bar Height:</strong> TMB contribution of each signature</li>
                <li><strong>Dot Size:</strong> Proportion of samples with the signature</li>
                <li><strong>Colors:</strong> Each signature has a specific color based on its biological etiology</li>
              </ul>
            </Col>
          </Row>
        </Card.Body>
      </Card>
    </div>
  );
}
