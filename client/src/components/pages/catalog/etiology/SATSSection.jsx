import React from 'react';
import { Row, Col, Card } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import Plot from 'react-plotly.js';
import { useSatsDataBySignatureQuery, useSatsEtiologyLookupQuery, useSatsExampleDataQuery } from './satsApiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

export default function SATSSection({ selectedSignature }) {

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
    
    // Keep signature names as they are - no mapping needed
    return signatureName;
  }

    // Function to get the correct signatureSetName for STS signatures
  function getSTSSignatureSetName(signatureName) {
    if (!signatureName) return null;
    
    if (signatureName.includes('SBS')) {
      // For SBS2_13, use the specific SATS signatureSetName
      if (signatureName === 'SBS2_13') {
        return 'SATS_TS_AACR_GENIE_GRCh37_SBS96';
      }
      return 'COSMIC_v3.4_Signatures_GRCh37_SBS96';
    } else if (signatureName.includes('DBS')) {
      // For DBS signatures, use the DBS78 signatureSetName
      return 'SATS_TS_AACR_GENIE_GRCh37_DBS78';
    } else if (signatureName.includes('ID')) {
      return 'COSMIC_v3.4_Signatures_GRCh37_ID83';
    }
    
    return null;
  }

  // Get mapped signature name for API calls
  const mappedSignatureName = mapSignatureName(signatureToUse);

  // Query the signature_etiology data to find signatureSetName for non-STS categories
  const { data: etiologyData, error: etiologyError, isLoading: etiologyLoading } = useSatsEtiologyLookupQuery(
    { limit: 1000 }, 
    { skip: !signatureToUse || category === 'STS' } // Skip for STS category
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

  // For STS category, use the hard-coded signatureSetName, otherwise use etiology lookup
  const signatureSetNameToUse = category === 'STS' 
    ? getSTSSignatureSetName(signatureToUse)
    : autoSignatureSetName;

  console.log('🎯 SATSSection signatureSetName resolution:', {
    category,
    signatureToUse,
    mappedSignatureName,
    autoSignatureSetName,
    stsSignatureSetName: getSTSSignatureSetName(signatureToUse),
    finalSignatureSetName: signatureSetNameToUse
  });

  // Only use example data for STS category (Signatures from Targeted Sequencing)
  const useExampleData = category === 'STS' && Boolean(signatureToUse);
  
  // Example data query (only run for STS category)
  const { 
    data: examplePlotConfig, 
    isFetching: fetchingExample, 
    error: exampleError 
  } = useSatsExampleDataQuery(signatureToUse, { skip: !useExampleData });

  // For non-STS categories, use real API data (skip when using example data)
  const { 
    data: plotConfig, 
    isFetching: fetchingPlot, 
    error: plotError 
  } = useSatsDataBySignatureQuery(
    { 
      signatureName: mappedSignatureName, // Use mapped signature name for API call
      signatureSetName: signatureSetNameToUse 
    },
    { skip: !signatureToUse || !signatureSetNameToUse || useExampleData }
  );

  // Choose which data to use based on category
  const finalPlotConfig = useExampleData ? examplePlotConfig : plotConfig;
  const finalFetching = useExampleData ? fetchingExample : fetchingPlot;
  const finalError = useExampleData ? exampleError : plotError;

  // Debug logging
  console.log('🎯 SATSSection Debug:', {
    category,
    signatureToUse,
    useExampleData,
    examplePlotConfig,
    plotConfig,
    finalPlotConfig,
    finalFetching,
    finalError
  });

  return (
    <div className="mt-4">
      <h5 className="separator">SATS - Signature Activity in Tumor Samples</h5>
      <p className="text-muted mb-3">
        Interactive visualization showing tumor mutational burden (TMB) and signature presence 
        across different cancer types for Signatures from Targeted Sequencing (STS). 
        This plot displays comprehensive signature activity data across all cancer types.
      </p>

      {/* Debug information for development */}
      {signatureToUse && (
        <div className="mb-3 p-2 bg-light rounded small">
          <strong>Analysis Info:</strong><br/>
          Selected Signature: {signatureToUse}<br/>
          Mapped for API: {mappedSignatureName}<br/>
          Signature Set: {signatureSetNameToUse || 'Determining...'}<br/>
          Category: {category}<br/>
          Data Source: {useExampleData ? 'Example Data' : 'API Data'}<br/>
        </div>
      )}

      <Card>
        <Card.Body style={{ padding: '1rem', display: 'flex', flexDirection: 'column', alignItems: 'center' }}>
          <LoadingOverlay active={finalFetching} />
          
          {finalError && (
            <div className="alert alert-danger">
              <h6>Error loading plot data</h6>
              <p>{finalError.message || 'An unknown error occurred'}</p>
            </div>
          )}

          {finalPlotConfig && finalPlotConfig.traces && finalPlotConfig.traces.length > 0 ? (
            <div style={{ display: 'flex', justifyContent: 'center', width: '100%', overflow: 'hidden' }}>
              <Plot
                data={JSON.parse(JSON.stringify(finalPlotConfig.traces))} // Deep copy to avoid read-only errors
                layout={JSON.parse(JSON.stringify(finalPlotConfig.layout))} // Deep copy layout too
                config={finalPlotConfig.config}
                style={{ width: '100%', height: '900px', maxWidth: '1400px' }}
                useResizeHandler={true}
                onInitialized={(figure, graphDiv) => console.log('✅ SATS Plot initialized successfully')}
                onError={(err) => console.error('❌ SATS Plot rendering error:', err)}
              />
            </div>
          ) : useExampleData && finalPlotConfig ? (
            <div className="text-center p-5">
              <h6>Plot Data Available But Not Rendering</h6>
              <p>finalPlotConfig exists but Plot component not showing</p>
              <pre style={{textAlign: 'left', fontSize: '10px'}}>
                {JSON.stringify({
                  hasTraces: !!finalPlotConfig.traces,
                  tracesLength: finalPlotConfig.traces?.length,
                  hasLayout: !!finalPlotConfig.layout,
                  hasConfig: !!finalPlotConfig.config
                }, null, 2)}
              </pre>
            </div>
          ) : signatureToUse && signatureSetNameToUse && !finalFetching && !finalError ? (
            <div className="text-center p-5 text-muted">
              <h6>No Data Available</h6>
              <p>No signature activity data found for {signatureToUse} in {signatureSetNameToUse}.</p>
              <p>This signature may not have activity data available in the current dataset.</p>
            </div>
          ) : signatureToUse && !signatureSetNameToUse ? (
            <div className="text-center p-5 text-muted">
              <h6>Analyzing Signature...</h6>
              <p>Determining signature set for {signatureToUse}...</p>
            </div>
          ) : !signatureToUse ? (
            <div className="text-center p-5 text-muted">
              <h6>Select a Signature</h6>
              <p>Choose a signature from the etiology section above to view the SATS plot.</p>
            </div>
          ) : null}
        </Card.Body>
      </Card>

      
    </div>
  );
}
