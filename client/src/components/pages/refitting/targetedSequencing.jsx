import React, { useState } from 'react';
import { Card, Row, Col, Table, Button, Form, Alert } from 'react-bootstrap';

export default function TargetedSequencing() {
  const [selectedMetric, setSelectedMetric] = useState('cosine');
  const [selectedAlgorithm, setSelectedAlgorithm] = useState('sigprofiler');

  // Mock data - this would typically come from props or API
  const mockResults = {
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

  const renderCosineResults = () => (
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
        {mockResults.cosineData.map((row, index) => (
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

  const renderBICResults = () => (
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
        {mockResults.bicData.map((row, index) => (
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

  const renderL2NormResults = () => (
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
        {mockResults.l2NormData.map((row, index) => (
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

  return (
    <div className="bg-white border rounded p-3">
      <h4>Targeted Sequencing Results</h4>
      
      <Alert variant="info">
        <strong>Analysis Complete:</strong> Your refitting analysis has been completed successfully. 
        Below are the comprehensive results showing the performance comparison between original and refitted signatures.
      </Alert>

      <Row className="mb-3">
        <Col md={6}>
          <Form.Group>
            <Form.Label><strong>Select Metric:</strong></Form.Label>
            <Form.Control
              as="select"
              value={selectedMetric}
              onChange={(e) => setSelectedMetric(e.target.value)}
            >
              <option value="cosine">Cosine Similarity</option>
              <option value="bic">BIC Comparison</option>
              <option value="l2norm">L2 Norm Analysis</option>
            </Form.Control>
          </Form.Group>
        </Col>
        <Col md={6}>
          <Form.Group>
            <Form.Label><strong>Algorithm Used:</strong></Form.Label>
            <Form.Control
              as="select"
              value={selectedAlgorithm}
              onChange={(e) => setSelectedAlgorithm(e.target.value)}
              disabled
            >
              <option value="sigprofiler">SigProfiler</option>
              <option value="deconstructsigs">deconstructSigs</option>
              <option value="bootstrap">Bootstrapping Method</option>
            </Form.Control>
          </Form.Group>
        </Col>
      </Row>

      <Card className="mb-4">
        <Card.Header>
          <h5 className="mb-0">
            {selectedMetric === 'cosine' && 'Cosine Similarity Analysis'}
            {selectedMetric === 'bic' && 'BIC Model Comparison'}
            {selectedMetric === 'l2norm' && 'L2 Norm Error Analysis'}
          </h5>
        </Card.Header>
        <Card.Body>
          {selectedMetric === 'cosine' && (
            <>
              <p>
                Cosine similarity measures how well the refitted signatures match the original signatures. 
                Values closer to 1.0 indicate better similarity.
              </p>
              {renderCosineResults()}
            </>
          )}
          
          {selectedMetric === 'bic' && (
            <>
              <p>
                BIC (Bayesian Information Criterion) comparison shows model performance. 
                Lower BIC scores and negative Δ AIC values indicate better model fit.
              </p>
              {renderBICResults()}
            </>
          )}
          
          {selectedMetric === 'l2norm' && (
            <>
              <p>
                L2 Norm analysis shows the reconstruction error for each sample. 
                Lower values and positive improvements indicate better refitting performance.
              </p>
              {renderL2NormResults()}
            </>
          )}
        </Card.Body>
      </Card>

      <Card>
        <Card.Header>
          <h5 className="mb-0">Summary & Recommendations</h5>
        </Card.Header>
        <Card.Body>
          <Row>
            <Col md={6}>
              <h6>Key Findings:</h6>
              <ul>
                <li>Average cosine similarity: <strong>0.96</strong></li>
                <li>BIC improvement: <strong>-44.4</strong></li>
                <li>Average L2 norm reduction: <strong>13.3%</strong></li>
                <li>Best performing algorithm: <strong>SigProfiler</strong></li>
              </ul>
            </Col>
            <Col md={6}>
              <h6>Recommendations:</h6>
              <ul>
                <li>SigProfiler shows the best overall performance</li>
                <li>All signatures show good refitting quality (similarity &gt; 0.94)</li>
                <li>Consider using the refitted signatures for downstream analysis</li>
                <li>Monitor SBS2 as it shows slightly lower performance</li>
              </ul>
            </Col>
          </Row>
          
          <div className="mt-3">
            <Button variant="primary" className="me-2">
              Download Results
            </Button>
            <Button variant="outline-secondary" className="me-2">
              Export Report
            </Button>
            <Button variant="outline-info">
              View Detailed Analysis
            </Button>
          </div>
        </Card.Body>
      </Card>
    </div>
  );
}
