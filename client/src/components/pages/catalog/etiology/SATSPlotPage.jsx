import React, { useState, useEffect } from 'react';
import { Form, Row, Col, Button, Card } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import Plot from 'react-plotly.js';
import { useSatsDataQuery } from './satsApiSlice';
import { useEtiologyOptionsQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Select from '../../../controls/select/selectHookForm';

export default function SATSPlotPage() {
  const [params, setParams] = useState(null);

  // Fetch etiology options for form controls
  const { data: options, isFetching: fetchingOptions } = useEtiologyOptionsQuery();

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
    <div className="p-3">
      <h2>Tissue Distribution</h2>
      <p className="lead">
        Interactive visualization showing tumor mutational burden (TMB) and signature presence 
        across different cancer types.
      </p>

      <Card className="mb-4">
        <Card.Header>
          <h5>Data Selection</h5>
        </Card.Header>
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
        <Card.Header>
          <h5>SATS Plot</h5>
        </Card.Header>
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
              style={{ width: '100%', height: '900px' }}
            />
          ) : params && !fetchingPlot && !plotError ? (
            <div className="text-center p-5 text-muted">
              <h5>No Data Available</h5>
              <p>No signature presence data found for the selected study and signature set.</p>
              <p>Try selecting a different study or signature set.</p>
            </div>
          ) : !params ? (
            <div className="text-center p-5 text-muted">
              <h5>Select Data to View Plot</h5>
              <p>Choose a study and signature set above to generate the SATS plot.</p>
            </div>
          ) : null}
        </Card.Body>
      </Card>

      
    </div>
  );
}
