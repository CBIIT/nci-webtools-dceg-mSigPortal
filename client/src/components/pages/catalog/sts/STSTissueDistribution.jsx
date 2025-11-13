import React from 'react';
import { Row, Col, Card } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import Plot from 'react-plotly.js';
import { useSatsDataBySignatureQuery, useSatsEtiologyLookupQuery, useSatsExampleDataQuery } from '../etiology/satsApiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

// Component to render a single SATS plot without header
function SATSPlot({ signatureName, signatureSetName, useExampleData, plotTitle, yAxisLabel }) {
  // For non-STS categories, use real API data (skip when using example data)
  const { 
    data: plotConfig, 
    isFetching: fetchingPlot, 
    error: plotError 
  } = useSatsDataBySignatureQuery(
    { 
      signatureName: signatureName,
      signatureSetName: signatureSetName 
    },
    { skip: !signatureName || !signatureSetName || useExampleData }
  );

  // Example data query (only run for STS category)
  const { 
    data: examplePlotConfig, 
    isFetching: fetchingExample, 
    error: exampleError 
  } = useSatsExampleDataQuery(signatureName, { skip: !useExampleData });

  // Choose which data to use based on category
  const finalPlotConfig = useExampleData ? examplePlotConfig : plotConfig;
  const finalFetching = useExampleData ? fetchingExample : fetchingPlot;
  const finalError = useExampleData ? exampleError : plotError;

  // Modify layout to include custom title and y-axis label
  const customLayout = finalPlotConfig?.layout ? {
    ...finalPlotConfig.layout,
    title: {
      text: plotTitle ? `<b>${plotTitle}</b>` : finalPlotConfig.layout.title?.text || finalPlotConfig.layout.title,
      font: {
        family: 'Times New Roman',
        size: 18
      },
      x: 0.5
    },
    yaxis: {
      ...finalPlotConfig.layout.yaxis,
      title: yAxisLabel || finalPlotConfig.layout.yaxis?.title
    }
  } : null;

  return (
    <Card className="mb-3">
      <Card.Body style={{ padding: '1rem', display: 'flex', flexDirection: 'column', alignItems: 'center' }}>
        <LoadingOverlay active={finalFetching} />
        
        {finalError && (
          <div className="alert alert-danger">
            <h6>Error loading plot data</h6>
            <p>{finalError.message || 'An unknown error occurred'}</p>
          </div>
        )}

        {finalPlotConfig && finalPlotConfig.traces && finalPlotConfig.traces.length > 0 && customLayout ? (
          <div style={{ display: 'flex', justifyContent: 'center', width: '100%', overflow: 'hidden' }}>
            <Plot
              data={JSON.parse(JSON.stringify(finalPlotConfig.traces))}
              layout={JSON.parse(JSON.stringify(customLayout))}
              config={finalPlotConfig.config}
              style={{ width: '100%', height: '900px' }}
              useResizeHandler={true}
            />
          </div>
        ) : signatureName && signatureSetName && !finalFetching && !finalError ? (
          <div className="text-center p-5 text-muted">
            <h6>No Data Available</h6>
            <p>No signature activity data found for {signatureName}.</p>
          </div>
        ) : !signatureName ? (
          <div className="text-center p-5 text-muted">
            <h6>Select a Signature</h6>
            <p>Choose a signature to view the SATS plot.</p>
          </div>
        ) : null}
      </Card.Body>
    </Card>
  );
}

export default function STSTissueDistribution({ selectedSignature }) {
  // Get current etiology category context
  const { category } = useSelector((state) => state.catalog.sts);

  // Determine signature set names
  function getSignatureSetName(signatureName) {
    if (!signatureName) return null;
    
    if (signatureName.includes('SBS')) {
      return 'SATS_TS_AACR_GENIE_GRCh37_SBS96';
    } else if (signatureName.includes('DBS')) {
      return 'SATS_TS_AACR_GENIE_GRCh37_DBS78';
    } else if (signatureName.includes('ID')) {
      return 'SATS_TS_AACR_GENIE_GRCh37_ID83';
    }
    
    return null;
  }

  // Determine which signatures to show
  const sbsSignature = selectedSignature?.includes('SBS') ? selectedSignature : 'SBS2_13';
  const dbsSignature = selectedSignature?.includes('DBS') ? selectedSignature : 'DBS2';

  const sbsSignatureSetName = getSignatureSetName(sbsSignature);
  const dbsSignatureSetName = getSignatureSetName(dbsSignature);

  return (
    <div className="mt-4">
      <h5 className="separator">Tissue Distribution</h5>
      
      {/* SBS Plot */}
      <SATSPlot 
        signatureName={sbsSignature}
        signatureSetName={sbsSignatureSetName}
        useExampleData={true}
        plotTitle="Repertoire of single base substitution (SBS) mutational signatures detected by targeted sequencing"
        yAxisLabel="Mutations per Mb"
      />

      {/* DBS Plot */}
      <SATSPlot 
        signatureName={dbsSignature}
        signatureSetName={dbsSignatureSetName}
        useExampleData={true}
        plotTitle="Repertoire of double base substitution (DBS) mutational signatures detected by targeted sequencing"
        yAxisLabel="Mutations per Mb"
      />
    </div>
  );
}
