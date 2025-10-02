import React from 'react';
import { Alert, ProgressBar, Table, Badge } from 'react-bootstrap';

export default function Status() {
  // Mock status data - this would typically come from props or API
  const mockStatus = {
    jobId: 'REF_20241002_001',
    status: 'IN_PROGRESS',
    submittedAt: '2024-10-02 14:30:00',
    estimatedCompletion: '2024-10-02 15:45:00',
    progress: 65,
    currentStep: 'Refitting signatures using SigProfiler',
    steps: [
      { name: 'Data validation', status: 'COMPLETED', timestamp: '2024-10-02 14:31:00' },
      { name: 'Signature extraction', status: 'COMPLETED', timestamp: '2024-10-02 14:35:00' },
      { name: 'Reference signature loading', status: 'COMPLETED', timestamp: '2024-10-02 14:40:00' },
      { name: 'Refitting with SigProfiler', status: 'IN_PROGRESS', timestamp: null },
      { name: 'Statistical evaluation', status: 'PENDING', timestamp: null },
      { name: 'Report generation', status: 'PENDING', timestamp: null }
    ]
  };

  const getStatusVariant = (status) => {
    switch (status) {
      case 'COMPLETED': return 'success';
      case 'IN_PROGRESS': return 'primary';
      case 'FAILED': return 'danger';
      case 'PENDING': return 'secondary';
      default: return 'info';
    }
  };

  const getStatusBadge = (status) => {
    switch (status) {
      case 'COMPLETED': return <Badge variant="success">Completed</Badge>;
      case 'IN_PROGRESS': return <Badge variant="primary">In Progress</Badge>;
      case 'FAILED': return <Badge variant="danger">Failed</Badge>;
      case 'PENDING': return <Badge variant="secondary">Pending</Badge>;
      default: return <Badge variant="info">Unknown</Badge>;
    }
  };

  return (
    <div className="bg-white border rounded p-3">
      <h4>Refitting Status</h4>
      
      <div className="mb-4">
        <div className="row">
          <div className="col-md-6">
            <strong>Job ID:</strong> {mockStatus.jobId}
          </div>
          <div className="col-md-6">
            <strong>Status:</strong> {getStatusBadge(mockStatus.status)}
          </div>
        </div>
        <div className="row mt-2">
          <div className="col-md-6">
            <strong>Submitted:</strong> {mockStatus.submittedAt}
          </div>
          <div className="col-md-6">
            <strong>Est. Completion:</strong> {mockStatus.estimatedCompletion}
          </div>
        </div>
      </div>

      {mockStatus.status === 'IN_PROGRESS' && (
        <div className="mb-4">
          <div className="d-flex justify-content-between align-items-center mb-2">
            <strong>Progress</strong>
            <span>{mockStatus.progress}%</span>
          </div>
          <ProgressBar 
            now={mockStatus.progress} 
            variant={getStatusVariant(mockStatus.status)}
            animated 
          />
          <small className="text-muted mt-2 d-block">
            Current step: {mockStatus.currentStep}
          </small>
        </div>
      )}

      {mockStatus.status === 'COMPLETED' && (
        <Alert variant="success">
          <strong>Analysis Complete!</strong> Your refitting analysis has finished successfully. 
          You can now view the results in the Targeted Sequencing tab.
        </Alert>
      )}

      {mockStatus.status === 'FAILED' && (
        <Alert variant="danger">
          <strong>Analysis Failed:</strong> There was an error processing your refitting analysis. 
          Please check your input data and try again, or contact support if the issue persists.
        </Alert>
      )}

      <h5>Processing Steps</h5>
      <Table striped bordered hover size="sm">
        <thead>
          <tr>
            <th>Step</th>
            <th>Status</th>
            <th>Timestamp</th>
          </tr>
        </thead>
        <tbody>
          {mockStatus.steps.map((step, index) => (
            <tr key={index}>
              <td>{step.name}</td>
              <td>{getStatusBadge(step.status)}</td>
              <td>{step.timestamp || '-'}</td>
            </tr>
          ))}
        </tbody>
      </Table>

      <div className="mt-3">
        <small className="text-muted">
          Note: This page will automatically refresh every minute while processing is in progress.
        </small>
      </div>
    </div>
  );
}
