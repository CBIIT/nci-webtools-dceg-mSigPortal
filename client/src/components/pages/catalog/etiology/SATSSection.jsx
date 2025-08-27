import React, { useState, useEffect } from 'react';
import { Form, Row, Col, Button, Card } from 'react-bootstrap';
import { useForm } from 'react-hook-form';
import { useSelector } from 'react-redux';
import Plot from 'react-plotly.js';
import { useSatsDataQuery, useSatsOptionsQuery, useSatsDataBySignatureQuery, useSatsEtiologyLookupQuery } from './satsApiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Select from '../../../controls/select/selectHookForm';

export default function SATSSection({ selectedSignature }) {
  const [params, setParams] = useState(null);

  // Get current etiology category context
  const { category, signature } = useSelector((state) => state.catalog.etiology);

  // Use the selected signature from props or state
  const signatureToUse = selectedSignature || signature;

  // Always log that SATSSection is rendering
  console.log('SATSSection rendering with:', {
    selectedSignature,
    signature,
    signatureToUse,
    category
  });

  // Helper function to map signature names for API calls
  function mapSignatureName(signatureName) {
    if (!signatureName) return signatureName;
    
    // For signatures like SBS2_13, map to SBS2
    const match = signatureName.match(/^(SBS\d+|DBS\d+|ID\d+)(_\d+)?$/);
    if (match) {
      return match[1]; // Return the base signature name without suffix
    }
    
    return signatureName;
  }

  // Get mapped signature name for API calls
  const mappedSignatureName = mapSignatureName(signatureToUse);

  // Fetch sample etiology data to get available studies and signature sets
  const { data: options, isFetching: fetchingOptions } = useSatsOptionsQuery();

  // Query the signature_etiology data to find signatureSetName for the selected signature
  const { data: etiologyData, error: etiologyError, isLoading: etiologyLoading } = useSatsEtiologyLookupQuery(
    { limit: 1000 }, 
    { skip: !signatureToUse }
  );

  // Log query status
  console.log('SATSSection etiology query status:', {
    signatureToUse,
    skip: !signatureToUse,
    isLoading: etiologyLoading,
    hasData: !!etiologyData,
    dataLength: etiologyData?.length || 0,
    error: etiologyError
  });

  // Find the signatureSetName for the selected signature from the etiology data
  const signatureInfo = etiologyData?.find(item => 
    item.signatureName === mappedSignatureName || 
    item.signature === mappedSignatureName ||
    item.signatureName === signatureToUse ||
    item.signature === signatureToUse
  );
  const autoSignatureSetName = signatureInfo?.signatureSetName;

  // Debug: Log the first few items to see the structure
  if (etiologyData && etiologyData.length > 0 && signatureToUse) {
    console.log('Etiology data sample:', etiologyData.slice(0, 3));
    console.log('Original signature:', signatureToUse);
    console.log('Mapped signature:', mappedSignatureName);
    console.log('Found signature info:', signatureInfo);
    
    // Log all unique signature names to see what's available
    const uniqueSignatureNames = [...new Set(etiologyData.map(item => item.signatureName || item.signature).filter(Boolean))].slice(0, 20);
    console.log('Available signature names (first 20):', uniqueSignatureNames);
    
    // Find all unique signatureSetNames for SBS signatures
    const sbsSignatureSetNames = [...new Set(
      etiologyData
        .filter(item => (item.signatureName || item.signature || '').includes('SBS'))
        .map(item => item.signatureSetName)
        .filter(Boolean)
    )];
    
    console.log('Available SBS signatureSetNames:', sbsSignatureSetNames);
    
    // Find all signatures that start with SBS2
    const sbs2Signatures = etiologyData.filter(item => 
      (item.signatureName || item.signature || '').startsWith('SBS2')
    );
    
    console.log('Available SBS2* signatures:', sbs2Signatures.map(item => ({
      signatureName: item.signatureName || item.signature,
      signatureSetName: item.signatureSetName
    })));
    
    // Show all unique signatureSetNames in the data
    const allSignatureSetNames = [...new Set(etiologyData.map(item => item.signatureSetName).filter(Boolean))];
    console.log('All unique signatureSetNames in etiology data:', allSignatureSetNames);
    
    // Look for signatures that contain 'SBS2'
    const sbs2Related = etiologyData.filter(item => 
      (item.signatureName && item.signatureName.includes('SBS2')) ||
      (item.signature && item.signature.includes('SBS2'))
    );
    console.log('SBS2-related signatures found:', sbs2Related.map(item => item.signatureName || item.signature));
  }

  // Automatically fetch SATS data when signature is selected
  const { 
    data: plotConfig, 
    isFetching: fetchingPlot, 
    error: plotError 
  } = useSatsDataBySignatureQuery(
    { 
      signatureName: mappedSignatureName, // Use mapped signature name for API call
      signatureSetName: autoSignatureSetName 
    },
    { skip: !signatureToUse || !autoSignatureSetName }
  );

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
        .sort((a, b) => {
          // Prioritize signature sets based on current category
          if (category === 'Cosmic' && a.includes('COSMIC')) return -1;
          if (category === 'Cosmic' && b.includes('COSMIC')) return 1;
          if (category === 'STS' && a.includes('STS')) return -1;
          if (category === 'STS' && b.includes('STS')) return 1;
          return a.localeCompare(b);
        })
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
      strategy: data.strategy?.value || 'WES', // Default to WES to match R analysis
      limit: 1000000 // Ensure we get all data
    };
    setParams(queryParams);
  }

  return (
    <div className="mt-4">
      <h5 className="separator">SATS - Signature Activity in Tumor Samples</h5>
      <p className="text-muted mb-3">
        Interactive visualization showing tumor mutational burden (TMB) and signature presence 
        across different cancer types. This plot reproduces the analysis from your R script, 
        showing both the TMB contribution per signature and the proportion of samples exhibiting each signature.
      </p>

      {/* Debug information for development */}
      {signatureToUse && (
        <div className="mb-3 p-2 bg-light rounded small">
          <strong>Debug Info:</strong><br/>
          Original Signature: {signatureToUse}<br/>
          Mapped Signature: {mappedSignatureName}<br/>
          Found SignatureSetName: {autoSignatureSetName || 'Not found'}<br/>
          Available Signatures: {etiologyData?.length || 0} total<br/>
          API Call: {autoSignatureSetName ? 
            `mutational_signature?signatureName=${mappedSignatureName}&signatureSetName=${autoSignatureSetName}` : 
            'Waiting for signatureSetName...'}
        </div>
      )}

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
          ) : signatureToUse && autoSignatureSetName && !fetchingPlot && !plotError ? (
            <div className="text-center p-5 text-muted">
              <h6>No Data Available</h6>
              <p>No signature activity data found for {signatureToUse} in {autoSignatureSetName}.</p>
              <p>This signature may not have activity data available in the current dataset.</p>
            </div>
          ) : signatureToUse && !autoSignatureSetName ? (
            <div className="text-center p-5 text-muted">
              <h6>Signature Set Not Found</h6>
              <p>Could not determine the signature set for {signatureToUse}.</p>
              <p>Use the form above to manually select a study and signature set.</p>
            </div>
          ) : params && !fetchingPlot && !plotError ? (
            <div className="text-center p-5 text-muted">
              <h6>No Data Available</h6>
              <p>No signature presence data found for the selected study and signature set.</p>
              <p>Try selecting a different study or signature set.</p>
            </div>
          ) : !params && !signatureToUse ? (
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
