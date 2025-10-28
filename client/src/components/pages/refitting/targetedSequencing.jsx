import React, { useState, useEffect } from 'react';
import { Card, Row, Col, Table, Button, Form, Alert } from 'react-bootstrap';

export default function TargetedSequencing({ jobId }) {
  const [selectedMetric, setSelectedMetric] = useState('h_est');
  const [selectedAlgorithm, setSelectedAlgorithm] = useState('sigprofiler');
  const [csvData, setCsvData] = useState([]);
  const [jobParams, setJobParams] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  // Load CSV data and job parameters when jobId changes
  useEffect(() => {
    if (jobId) {
      loadJobData(jobId);
    }
  }, [jobId]);

  const loadJobData = async (jobId) => {
    setLoading(true);
    setError(null);
    try {
      console.log(`Loading job data for: ${jobId}`);
      
      // Load job parameters from params.json
      const paramsResponse = await fetch(`/mutational-signatures/api/data/input/${jobId}/params.json`);
      if (paramsResponse.ok) {
        const paramsText = await paramsResponse.text();
        const params = JSON.parse(paramsText);
        console.log('Loaded job parameters:', params);
        setJobParams(params);
      } else {
        console.warn('Could not load job parameters, using defaults');
        setJobParams({ signatureType: 'SBS' }); // Default fallback
      }
      
      // Load CSV data
      const csvResponse = await fetch(`/mutational-signatures/api/data/output/${jobId}/H_Burden_est.csv`);
      if (!csvResponse.ok) {
        throw new Error(`Failed to load CSV data: ${csvResponse.status} ${csvResponse.statusText}`);
      }
      const csvText = await csvResponse.text();
      console.log('Raw CSV text:', csvText.substring(0, 500)); // Log first 500 chars
      const parsedData = parseCsv(csvText);
      console.log('Parsed CSV data:', parsedData.slice(0, 5)); // Log first 5 rows
      setCsvData(parsedData);
    } catch (err) {
      setError(err.message);
      console.error('Error loading job data:', err);
    } finally {
      setLoading(false);
    }
  };

  const parseCsv = (csvText) => {
    const lines = csvText.trim().split('\n');
    const headers = lines[0].split(',').map(h => h.trim().replace(/"/g, ''));
    
    return lines.slice(1).map(line => {
      // Handle quoted CSV values properly
      const values = line.split(',').map(v => v.trim().replace(/"/g, ''));
      const row = {};
      headers.forEach((header, index) => {
        row[header] = values[index] || '';
      });
      return row;
    });
  };

  // Get signature type from job parameters
  const submittedSignatureType = jobParams?.signatureType || 'SBS';

  const getCurrentResults = () => {
    // If we have real CSV data and a jobId, use it
    if (jobId && csvData.length > 0) {
      return {
        h_estData: csvData.map(row => ({
          sample_id: row.SAMPLE_ID || row.sample_id || '',
          signature: row.Signature || row.signature || '',
          activity: parseFloat(row.Activity || row.activity || 0),
          burden: parseFloat(row.Burden || row.burden || 0)
        })),
        // For now, only h_estData is available from the CSV
        // Other data types would come from additional CSV files or API endpoints
        burden_estData: [],
        cosineData: [],
        bicData: [],
        l2NormData: []
      };
    }
    
    // Return empty structure when no data is available
    return {
      h_estData: [],
      burden_estData: [],
      cosineData: [],
      bicData: [],
      l2NormData: []
    };
  };

  const renderHEstResults = () => {
    const results = getCurrentResults();
    
    // Get unique samples and signatures
    const samples = [...new Set(results.h_estData.map(row => row.sample_id))];
    const signatures = [...new Set(results.h_estData.map(row => row.signature))].sort();
    
    // Create a lookup object for easy data access
    const dataLookup = {};
    results.h_estData.forEach(row => {
      if (!dataLookup[row.sample_id]) {
        dataLookup[row.sample_id] = {};
      }
      dataLookup[row.sample_id][row.signature] = {
        activity: row.activity,
        burden: row.burden
      };
    });
    
    return (
      <div>
        <h6>Activity Values</h6>
        <Table striped bordered hover responsive className="mb-4">
          <thead>
            <tr>
              <th>H_est</th>
              {signatures.map(sig => (
                <th key={sig}>{sig}</th>
              ))}
            </tr>
          </thead>
          <tbody>
            {samples.map(sample => (
              <tr key={sample}>
                <td><strong>{sample}</strong></td>
                {signatures.map(sig => (
                  <td key={sig}>
                    {dataLookup[sample]?.[sig]?.activity ?? 'N/A'}
                  </td>
                ))}
              </tr>
            ))}
          </tbody>
        </Table>
        
        <h6>Burden Values</h6>
        <Table striped bordered hover responsive>
          <thead>
            <tr>
              <th>Burden_est</th>
              {signatures.map(sig => (
                <th key={sig}>{sig}</th>
              ))}
            </tr>
          </thead>
          <tbody>
            {samples.map(sample => (
              <tr key={sample}>
                <td><strong>{sample}</strong></td>
                {signatures.map(sig => (
                  <td key={sig}>
                    {dataLookup[sample]?.[sig]?.burden ?? 'N/A'}
                  </td>
                ))}
              </tr>
            ))}
          </tbody>
        </Table>
      </div>
    );
  };

  const renderBurdenEstResults = () => {
    const results = getCurrentResults();
    
    if (!results.burden_estData || results.burden_estData.length === 0) {
      return (
        <Alert variant="info">
          <strong>No Burden Estimation Data:</strong> This data is not available in the current H_Burden_est.csv file. 
          Additional analysis outputs would be needed to display burden estimation results.
        </Alert>
      );
    }
    
    return (
      <Table striped bordered hover>
        <thead>
          <tr>
            <th>Sample</th>
            <th>Total Mutations</th>
            <th>Estimated Burden</th>
            <th>Accuracy</th>
          </tr>
        </thead>
        <tbody>
          {results.burden_estData.map((row, index) => (
            <tr key={index}>
              <td><strong>{row.sample}</strong></td>
              <td>{row.total_mutations.toLocaleString()}</td>
              <td>{row.estimated_burden.toLocaleString()}</td>
              <td>
                <span className="text-success">
                  {row.accuracy}
                </span>
              </td>
            </tr>
          ))}
        </tbody>
      </Table>
    );
  };

  const renderCosineResults = () => {
    const results = getCurrentResults();
    
    if (!results.cosineData || results.cosineData.length === 0) {
      return (
        <Alert variant="info">
          <strong>No Cosine Similarity Data:</strong> This analysis is not available in the current H_Burden_est.csv file. 
          Additional analysis outputs would be needed to display cosine similarity results.
        </Alert>
      );
    }
    
    return (
      <Table striped bordered hover>
        <thead>
          <tr>
            <th>Signature</th>
            <th>Original Weight</th>
            <th>Refitted Weight</th>
            <th>Cosine Similarity</th>
          </tr>
        </thead>
        <tbody>
          {results.cosineData.map((row, index) => (
            <tr key={index}>
              <td><strong>{row.signature}</strong></td>
              <td>{row.original.toFixed(3)}</td>
              <td>{row.refitted.toFixed(3)}</td>
              <td>
                <span className={row.similarity >= 0.95 ? 'text-success' : row.similarity >= 0.90 ? 'text-warning' : 'text-danger'}>
                  {row.similarity.toFixed(3)}
                </span>
              </td>
            </tr>
          ))}
        </tbody>
      </Table>
    );
  };

  const renderBICResults = () => {
    const results = getCurrentResults();
    
    if (!results.bicData || results.bicData.length === 0) {
      return (
        <Alert variant="info">
          <strong>No BIC Analysis Data:</strong> This analysis is not available in the current H_Burden_est.csv file. 
          Additional analysis outputs would be needed to display BIC model comparison results.
        </Alert>
      );
    }
    
    return (
      <Table striped bordered hover>
        <thead>
          <tr>
            <th>Model</th>
            <th>Signatures</th>
            <th>BIC Score</th>
            <th>Δ AIC</th>
          </tr>
        </thead>
        <tbody>
          {results.bicData.map((row, index) => (
            <tr key={index}>
              <td><strong>{row.model}</strong></td>
              <td>{row.signatures}</td>
              <td>{row.bic.toFixed(1)}</td>
              <td>
                <span className={row.deltaAIC < 0 ? 'text-success' : 'text-muted'}>
                  {row.deltaAIC.toFixed(1)}
                </span>
              </td>
            </tr>
          ))}
        </tbody>
      </Table>
    );
  };

  const renderL2NormResults = () => {
    const results = getCurrentResults();
    
    if (!results.l2NormData || results.l2NormData.length === 0) {
      return (
        <Alert variant="info">
          <strong>No L2 Norm Analysis Data:</strong> This analysis is not available in the current H_Burden_est.csv file. 
          Additional analysis outputs would be needed to display L2 norm error analysis results.
        </Alert>
      );
    }
    
    return (
      <Table striped bordered hover>
        <thead>
          <tr>
            <th>Sample</th>
            <th>Original L2 Norm</th>
            <th>Refitted L2 Norm</th>
            <th>Improvement</th>
          </tr>
        </thead>
        <tbody>
          {results.l2NormData.map((row, index) => (
            <tr key={index}>
              <td><strong>{row.sample}</strong></td>
              <td>{row.original.toFixed(3)}</td>
              <td>{row.refitted.toFixed(3)}</td>
              <td>
                <span className="text-success">
                  {row.improvement}
                </span>
              </td>
            </tr>
          ))}
        </tbody>
      </Table>
    );
  };

  return (
    <div className="bg-white border rounded p-3">
      <h4>Targeted Sequencing Results</h4>
      
      {jobParams && (
        <div className="mb-3">
          <small className="text-muted">
            <strong>Job:</strong> {jobParams.jobName || jobParams.jobId} | 
            <strong> Signature Type:</strong> {jobParams.signatureType} | 
            <strong> Genome:</strong> {jobParams.genome}
          </small>
        </div>
      )}
      
      {!jobId && (
        <Alert variant="warning">
          <strong>No Job Selected:</strong> Please select a job from the Status tab to view targeted sequencing results.
        </Alert>
      )}
      
      {jobId && loading && (
        <Alert variant="info">
          <strong>Loading:</strong> Loading analysis results for job {jobId}...
        </Alert>
      )}
      
      {jobId && error && (
        <Alert variant="danger">
          <strong>Error:</strong> {error}
        </Alert>
      )}
      
      {jobId && !loading && !error && (
        <>
          <Alert 
            variant="info" 
            style={{ 
              backgroundColor: '#689f39', 
              borderColor: '#689f39', 
              color: 'white' 
            }}
          >
            <strong>Analysis Complete:</strong> Your refitting analysis has been completed successfully.         
          </Alert>

      <Card className="mb-4">
        <Card.Header>
          <h5 className="mb-0">
            {submittedSignatureType} - {selectedMetric === 'h_est' && 'Refitting Results (H_Burden_est.csv)'}
            {selectedMetric === 'burden_est' && 'Burden_est Table - Mutation Burden Estimation'}
            {selectedMetric === 'cosine' && 'Cosine Similarity Analysis'}
            {selectedMetric === 'bic' && 'BIC Model Comparison'}
            {selectedMetric === 'l2norm' && 'L2 Norm Error Analysis'}
          </h5>
        </Card.Header>
        <Card.Body>
          {selectedMetric === 'h_est' && (
            <>
              {renderHEstResults()}
            </>
          )}

          {selectedMetric === 'burden_est' && (
            <>
              <p>
                The Burden_est table shows the estimated mutational burden for each sample compared to the total observed mutations. 
                Accuracy indicates how well the refitting process captured the original mutational burden.
              </p>
              {renderBurdenEstResults()}
            </>
          )}
          
          {selectedMetric === 'cosine' && (
            <>
              <p>
                Cosine similarity measures how well the refitted {submittedSignatureType} signatures match the original signatures. 
                Values closer to 1.0 indicate better similarity.
              </p>
              {renderCosineResults()}
            </>
          )}
          
          {selectedMetric === 'bic' && (
            <>
              <p>
                BIC (Bayesian Information Criterion) comparison shows model performance for {submittedSignatureType} signatures. 
                Lower BIC scores and negative Δ AIC values indicate better model fit.
              </p>
              {renderBICResults()}
            </>
          )}
          
          {selectedMetric === 'l2norm' && (
            <>
              <p>
                L2 Norm analysis shows the reconstruction error for each sample using {submittedSignatureType} signatures. 
                Lower values and positive improvements indicate better refitting performance.
              </p>
              {renderL2NormResults()}
            </>
          )}
        </Card.Body>
      </Card>
        </>
      )}

      
    </div>
  );
}
