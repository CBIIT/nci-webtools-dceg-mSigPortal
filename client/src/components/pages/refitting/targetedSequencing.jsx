import React, { useState } from 'react';
import { Card, Row, Col, Table, Button, Form, Alert, Nav, Tab } from 'react-bootstrap';

export default function TargetedSequencing() {
  const [selectedMetric, setSelectedMetric] = useState('h_est');
  const [selectedAlgorithm, setSelectedAlgorithm] = useState('sigprofiler');
  const [activeTab, setActiveTab] = useState('sbs');

  // Mock signature type from submitted job - this would come from props/context in real implementation
  const submittedSignatureType = 'SBS'; // This would be passed from the form submission

  // Mock data for SBS results - based on actual R function output
  const sbsResults = {
    h_estData: [
      { sample_id: 'GENIE-DFCI-050984-218969', signature: 'SBS1', activity: 0, burden: 0 },
      { sample_id: 'GENIE-DFCI-050984-218969', signature: 'SBS10a', activity: 0, burden: 0 },
      { sample_id: 'GENIE-DFCI-050984-218969', signature: 'SBS10b', activity: 0, burden: 0 },
      { sample_id: 'GENIE-DFCI-050984-218969', signature: 'SBS40a', activity: 0.003, burden: 0 },
      { sample_id: 'GENIE-DFCI-050984-218969', signature: 'SBS5', activity: 28.034, burden: 2 },
      { sample_id: 'GENIE-DFCI-050984-218969', signature: 'SBS6', activity: 0.001, burden: 0 },
      { sample_id: 'GENIE-DFCI-109295-436549', signature: 'SBS1', activity: 0, burden: 0 },
      { sample_id: 'GENIE-DFCI-109295-436549', signature: 'SBS10b', activity: 0, burden: 0 },
      { sample_id: 'GENIE-DFCI-109295-436549', signature: 'SBS5', activity: 32.735, burden: 2.34 },
      { sample_id: 'GENIE-DFCI-109295-436549', signature: 'SBS7a', activity: 37.782, burden: 2.65 },
      { sample_id: 'GENIE-DFCI-109295-436549', signature: 'SBS7b', activity: 0.16, burden: 0.01 }
    ],
    burden_estData: [
      { 
        sample: 'GENIE-DFCI-050984-218969', 
        total_mutations: 1847, estimated_burden: 1623, accuracy: '87.9%' 
      },
      { 
        sample: 'GENIE-DFCI-109295-436549', 
        total_mutations: 2156, estimated_burden: 1934, accuracy: '89.7%' 
      }
    ],
    cosineData: [
      { signature: 'SBS1', original: 0.95, refitted: 0.92, similarity: 0.97 },
      { signature: 'SBS2', original: 0.88, refitted: 0.85, similarity: 0.94 },
      { signature: 'SBS3', original: 0.76, refitted: 0.79, similarity: 0.96 },
      { signature: 'SBS13', original: 0.82, refitted: 0.80, similarity: 0.98 }
    ],
    bicData: [
      { model: 'Original', signatures: 4, bic: -2847.3, deltaAIC: 0.0 },
      { model: 'Refitted (SigProfiler)', signatures: 4, bic: -2891.7, deltaAIC: -44.4 },
      { model: 'Refitted (deconstructSigs)', signatures: 3, bic: -2856.2, deltaAIC: -8.9 }
    ],
    l2NormData: [
      { sample: 'Sample_001', original: 0.23, refitted: 0.19, improvement: '17.4%' },
      { sample: 'Sample_002', original: 0.31, refitted: 0.28, improvement: '9.7%' },
      { sample: 'Sample_003', original: 0.18, refitted: 0.15, improvement: '16.7%' },
      { sample: 'Sample_004', original: 0.42, refitted: 0.38, improvement: '9.5%' }
    ]
  };

  // Mock data for DBS results - following same format as SBS
  const dbsResults = {
    h_estData: [
      { sample_id: 'GENIE-DFCI-050984-218969', signature: 'DBS1', activity: 0, burden: 0 },
      { sample_id: 'GENIE-DFCI-050984-218969', signature: 'DBS2', activity: 12.456, burden: 1.2 },
      { sample_id: 'GENIE-DFCI-050984-218969', signature: 'DBS6', activity: 0.002, burden: 0 },
      { sample_id: 'GENIE-DFCI-050984-218969', signature: 'DBS11', activity: 8.923, burden: 0.8 },
      { sample_id: 'GENIE-DFCI-109295-436549', signature: 'DBS1', activity: 15.234, burden: 1.5 },
      { sample_id: 'GENIE-DFCI-109295-436549', signature: 'DBS2', activity: 23.567, burden: 2.1 },
      { sample_id: 'GENIE-DFCI-109295-436549', signature: 'DBS6', activity: 0.045, burden: 0 },
      { sample_id: 'GENIE-DFCI-109295-436549', signature: 'DBS11', activity: 18.789, burden: 1.7 }
    ],
    burden_estData: [
      { 
        sample: 'GENIE-DFCI-050984-218969', 
        total_mutations: 342, estimated_burden: 298, accuracy: '87.1%' 
      },
      { 
        sample: 'GENIE-DFCI-109295-436549', 
        total_mutations: 458, estimated_burden: 412, accuracy: '90.0%' 
      }
    ],
    cosineData: [
      { signature: 'DBS1', original: 0.89, refitted: 0.87, similarity: 0.95 },
      { signature: 'DBS2', original: 0.93, refitted: 0.91, similarity: 0.97 },
      { signature: 'DBS6', original: 0.71, refitted: 0.74, similarity: 0.94 },
      { signature: 'DBS11', original: 0.86, refitted: 0.84, similarity: 0.96 }
    ],
    bicData: [
      { model: 'Original', signatures: 4, bic: -1923.7, deltaAIC: 0.0 },
      { model: 'Refitted (SigProfiler)', signatures: 4, bic: -1967.2, deltaAIC: -43.5 },
      { model: 'Refitted (deconstructSigs)', signatures: 3, bic: -1931.4, deltaAIC: -7.7 }
    ],
    l2NormData: [
      { sample: 'Sample_001', original: 0.27, refitted: 0.23, improvement: '14.8%' },
      { sample: 'Sample_002', original: 0.34, refitted: 0.31, improvement: '8.8%' },
      { sample: 'Sample_003', original: 0.21, refitted: 0.18, improvement: '14.3%' },
      { sample: 'Sample_004', original: 0.39, refitted: 0.35, improvement: '10.3%' }
    ]
  };

  const getCurrentResults = () => {
    return activeTab === 'sbs' ? sbsResults : dbsResults;
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
      
      <Alert 
        variant="info" 
        style={{ 
          backgroundColor: '#689f39', 
          borderColor: '#689f39', 
          color: 'white' 
        }}
      >
        <strong>Analysis Complete:</strong> Your refitting analysis has been completed successfully. 
        Below are the comprehensive results showing the performance comparison between original and refitted signatures.
      </Alert>

      {/* SBS/DBS Tabs */}
      <Nav variant="tabs" className="mb-3">
        <Nav.Item>
          <Nav.Link 
            active={activeTab === 'sbs'}
            disabled={submittedSignatureType !== 'SBS'}
            onClick={() => setActiveTab('sbs')}
            style={{ 
              cursor: submittedSignatureType !== 'SBS' ? 'not-allowed' : 'pointer',
              opacity: submittedSignatureType !== 'SBS' ? 0.5 : 1,
              backgroundColor: activeTab === 'sbs' ? '#689f39 !important' : '#e8f5e8',
              color: activeTab === 'sbs' ? 'white !important' : '#689f39',
              borderColor: '#689f39 !important',
              fontWeight: activeTab === 'sbs' ? 'bold' : 'normal'
            }}
          >
            SBS Results
          </Nav.Link>
        </Nav.Item>
        <Nav.Item>
          <Nav.Link 
            active={activeTab === 'dbs'}
            disabled={submittedSignatureType !== 'DBS'}
            onClick={() => setActiveTab('dbs')}
            style={{ 
              cursor: submittedSignatureType !== 'DBS' ? 'not-allowed' : 'pointer',
              opacity: submittedSignatureType !== 'DBS' ? 0.5 : 1,
              backgroundColor: activeTab === 'dbs' ? '#689f39 !important' : '#e8f5e8',
              color: activeTab === 'dbs' ? 'white !important' : '#689f39',
              borderColor: '#689f39 !important',
              fontWeight: activeTab === 'dbs' ? 'bold' : 'normal'
            }}
          >
            DBS Results
          </Nav.Link>
        </Nav.Item>
      </Nav>

      

      <Card className="mb-4">
        <Card.Header>
          <h5 className="mb-0">
            {activeTab.toUpperCase()} - {selectedMetric === 'h_est' && 'SBS Refitting Results (H_Burden_est.csv)'}
            {selectedMetric === 'burden_est' && 'Burden_est Table - Mutation Burden Estimation'}
            {selectedMetric === 'cosine' && 'Cosine Similarity Analysis'}
            {selectedMetric === 'bic' && 'BIC Model Comparison'}
            {selectedMetric === 'l2norm' && 'L2 Norm Error Analysis'}
          </h5>
        </Card.Header>
        <Card.Body>
          {selectedMetric === 'h_est' && (
            <>
              <p>
                <strong>SBS Refitting Results:</strong> This table shows the output from the <code>run_sbs_refitting()</code> R function. 
                Each row represents a sample-signature combination with Activity (signature contribution strength) and Burden (estimated mutation count).
                Results are saved as "H_Burden_est.csv" and show significant signatures with non-zero activity values.
              </p>
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
                Cosine similarity measures how well the refitted {activeTab.toUpperCase()} signatures match the original signatures. 
                Values closer to 1.0 indicate better similarity.
              </p>
              {renderCosineResults()}
            </>
          )}
          
          {selectedMetric === 'bic' && (
            <>
              <p>
                BIC (Bayesian Information Criterion) comparison shows model performance for {activeTab.toUpperCase()} signatures. 
                Lower BIC scores and negative Δ AIC values indicate better model fit.
              </p>
              {renderBICResults()}
            </>
          )}
          
          {selectedMetric === 'l2norm' && (
            <>
              <p>
                L2 Norm analysis shows the reconstruction error for each sample using {activeTab.toUpperCase()} signatures. 
                Lower values and positive improvements indicate better refitting performance.
              </p>
              {renderL2NormResults()}
            </>
          )}
        </Card.Body>
      </Card>

      
    </div>
  );
}
