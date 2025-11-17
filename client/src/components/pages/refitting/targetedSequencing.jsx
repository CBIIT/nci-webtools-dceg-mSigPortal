import React, { useState, useEffect } from 'react';
import { Card, Row, Col, Button, Form, Alert } from 'react-bootstrap';
import { parseCSV } from '../../../services/utils';
import Table from '../../controls/table/table2';

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
      const parsedData = await parseCSV(csvText);
      console.log('Parsed CSV data:', parsedData.slice(0, 5)); // Log first 5 rows
      setCsvData(parsedData);
    } catch (err) {
      setError(err.message);
      console.error('Error loading job data:', err);
    } finally {
      setLoading(false);
    }
  };

  // Get signature type from job parameters
  const submittedSignatureType = jobParams?.signatureType || 'SBS';

  const getCurrentResults = () => {
    // If we have real CSV data and a jobId, use it
    if (jobId && csvData.length > 0) {
      return {
        h_estData: csvData
          .filter(row => row.SAMPLE_ID ) // Filter out empty rows
          .map(row => ({
            sample_id: row.SAMPLE_ID || '',
            signature: row.Signature || '',
            activity: parseFloat(row.Activity || 0),
            burden: parseFloat(row.Burden || 0)
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
    
    // Prepare data for Activity table
    const activityData = samples.map(sample => {
      const row = { sample_id: sample };
      signatures.forEach(sig => {
        row[sig] = dataLookup[sample]?.[sig]?.activity ?? 'N/A';
      });
      return row;
    });
    
    // Prepare data for Burden table
    const burdenData = samples.map(sample => {
      const row = { sample_id: sample };
      signatures.forEach(sig => {
        row[sig] = dataLookup[sample]?.[sig]?.burden ?? 'N/A';
      });
      return row;
    });
    
    // Define columns for Activity table
    const activityColumns = [
      {
        accessor: 'sample_id',
        Header: 'H_est',
        Cell: ({ value }) => <strong>{value}</strong>,
      },
      ...signatures.map(sig => ({
        accessor: sig,
        Header: sig,
      }))
    ];
    
    // Define columns for Burden table
    const burdenColumns = [
      {
        accessor: 'sample_id',
        Header: 'Burden_est',
        Cell: ({ value }) => <strong>{value}</strong>,
      },
      ...signatures.map(sig => ({
        accessor: sig,
        Header: sig,
      }))
    ];
    
    return (
      <div>
        <h3 className='h6-title'>Activity Values</h3>
        <Table 
          columns={activityColumns}
          data={activityData}
          striped
          bordered
        />

        <h3 className="mt-4 h6-title">Burden Values</h3>
        <Table 
          columns={burdenColumns}
          data={burdenData}
          striped
          bordered
        />
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
    
    const columns = [
      {
        accessor: 'sample',
        Header: 'Sample',
        Cell: ({ value }) => <strong>{value}</strong>,
      },
      {
        accessor: 'total_mutations',
        Header: 'Total Mutations',
        Cell: ({ value }) => value.toLocaleString(),
      },
      {
        accessor: 'estimated_burden',
        Header: 'Estimated Burden',
        Cell: ({ value }) => value.toLocaleString(),
      },
      {
        accessor: 'accuracy',
        Header: 'Accuracy',
        Cell: ({ value }) => (
          <span className="text-success">{value}</span>
        ),
      },
    ];
    
    return (
      <Table 
        columns={columns}
        data={results.burden_estData}
        striped
        bordered
      />
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
    
    const columns = [
      {
        accessor: 'signature',
        Header: 'Signature',
        Cell: ({ value }) => <strong>{value}</strong>,
      },
      {
        accessor: 'original',
        Header: 'Original Weight',
        Cell: ({ value }) => value.toFixed(3),
      },
      {
        accessor: 'refitted',
        Header: 'Refitted Weight',
        Cell: ({ value }) => value.toFixed(3),
      },
      {
        accessor: 'similarity',
        Header: 'Cosine Similarity',
        Cell: ({ value }) => (
          <span className={value >= 0.95 ? 'text-success' : value >= 0.90 ? 'text-warning' : 'text-danger'}>
            {value.toFixed(3)}
          </span>
        ),
      },
    ];
    
    return (
      <Table 
        columns={columns}
        data={results.cosineData}
        striped
        bordered
      />
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
    
    const columns = [
      {
        accessor: 'model',
        Header: 'Model',
        Cell: ({ value }) => <strong>{value}</strong>,
      },
      {
        accessor: 'signatures',
        Header: 'Signatures',
      },
      {
        accessor: 'bic',
        Header: 'BIC Score',
        Cell: ({ value }) => value.toFixed(1),
      },
      {
        accessor: 'deltaAIC',
        Header: 'Δ AIC',
        Cell: ({ value }) => (
          <span className={value < 0 ? 'text-success' : 'text-muted'}>
            {value.toFixed(1)}
          </span>
        ),
      },
    ];
    
    return (
      <Table 
        columns={columns}
        data={results.bicData}
        striped
        bordered
      />
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
    
    const columns = [
      {
        accessor: 'sample',
        Header: 'Sample',
        Cell: ({ value }) => <strong>{value}</strong>,
      },
      {
        accessor: 'original',
        Header: 'Original L2 Norm',
        Cell: ({ value }) => value.toFixed(3),
      },
      {
        accessor: 'refitted',
        Header: 'Refitted L2 Norm',
        Cell: ({ value }) => value.toFixed(3),
      },
      {
        accessor: 'improvement',
        Header: 'Improvement',
        Cell: ({ value }) => (
          <span className="text-success">{value}</span>
        ),
      },
    ];
    
    return (
      <Table 
        columns={columns}
        data={results.l2NormData}
        striped
        bordered
      />
    );
  };

  return (
    <div className="bg-white border rounded p-3">
      <h1 className='h4-title'>Targeted Sequencing Results</h1>
      
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
              backgroundColor: '#406024', 
              borderColor: '#406024', 
              color: 'white' 
            }}
          >
            <strong>Analysis Complete:</strong> Your refitting analysis has been completed successfully.         
          </Alert>

      <Card className="mb-4">
        <Card.Header>
          <h2 className="h5-title mb-0">
            {submittedSignatureType} - {selectedMetric === 'h_est' && 'Refitting Results (H_Burden_est.csv)'}
            {selectedMetric === 'burden_est' && 'Burden_est Table - Mutation Burden Estimation'}
            {selectedMetric === 'cosine' && 'Cosine Similarity Analysis'}
            {selectedMetric === 'bic' && 'BIC Model Comparison'}
            {selectedMetric === 'l2norm' && 'L2 Norm Error Analysis'}
          </h2>
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
