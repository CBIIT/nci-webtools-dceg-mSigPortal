import { useState, useMemo, useEffect } from 'react';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { Button, Popover, OverlayTrigger } from 'react-bootstrap';
import { Container } from 'react-bootstrap';
import { Link } from 'react-router-dom';
import Table from '../../controls/table/table2';
import moment from 'moment';
import momentDurationFormatSetup from 'moment-duration-format';
momentDurationFormatSetup(moment);

export default function Status() {
  const [jobs, setJobs] = useState([]);
  const [data, setData] = useState([]);
  const [isFetching, setIsFetching] = useState(false);

  // get jobs from local storage
  useEffect(() => {
    const localJobs = JSON.parse(localStorage.getItem('refitting-jobs') || '[]');
    console.log('Loading refitting jobs from localStorage:', localJobs);
    setJobs(localJobs);
    setData(localJobs); // Use stored jobs as data for now
  }, []);

  // update jobs in local storage
  useEffect(() => {
    if (jobs && jobs.length) {
      localStorage.setItem('refitting-jobs', JSON.stringify(jobs));
      setData(jobs);
      console.log('Updated refitting jobs in localStorage:', jobs);
    }
  }, [jobs]);

  const refetch = () => {
    // Simulate refresh by reloading from localStorage
    const localJobs = JSON.parse(localStorage.getItem('refitting-jobs') || '[]');
    console.log('Refetching refitting jobs:', localJobs);
    setJobs(localJobs);
    setData(localJobs);
  };

  function removeJob(id) {
    const remaining = jobs.filter((e) => e.id !== id);
    setJobs(remaining);
  }

  const columns = useMemo(() => [
    {
      accessor: 'jobName',
      Header: 'Name',
      Cell: ({ row, value }) => (
        <Link to={`/refitting/${row.original.id}`} style={{ textDecoration: 'none', color: '#007bff' }}>
          {value || 'Refitting Job'}
        </Link>
      ),
    },
    {
      accessor: 'status',
      Header: 'Status',
      Cell: ({ value }) => {
        if (value === 'SUBMITTED') return 'Submitted';
        else if (value === 'IN_PROGRESS') return 'In Progress';
        else if (value === 'FAILED') return 'Failed';
        else if (value === 'COMPLETED') return 'Completed';
        else return value || 'Submitted';
      },
    },
    {
      accessor: 'submittedAt',
      Header: 'Submitted Time',
      Cell: ({ value }) => {
        if (!value) return 'Just now';
        const time = moment(value);
        return `${time.format('LLL')} (${time.fromNow()})`;
      },
    },
    {
      accessor: 'id',
      Header: 'Job ID',
      Cell: ({ value }) => value || 'N/A',
    },
    {
      id: 'Download',
      accessor: 'status',
      Header: 'Download',
      Cell: ({ row, value }) => (
        <Button
          variant="link"
          href={`api/refitting/download/${row.original.id}/results.zip`}
          disabled={!['COMPLETED'].includes(value)}
        >
          Download
        </Button>
      ),
    },
    {
      Header: 'Remove',
      Cell: (e) => (
        <OverlayTrigger
          trigger="click"
          placement="left"
          rootClose
          overlay={
            <Popover>
              <Popover.Title as="h3">Confirm Remove</Popover.Title>
              <Popover.Content className="d-flex">
                <Button
                  className="w-100 mx-auto"
                  variant="danger"
                  onClick={() => removeJob(e.row.original.id)}
                >
                  Confirm
                </Button>
              </Popover.Content>
            </Popover>
          }
        >
          <Button variant="danger">Remove</Button>
        </OverlayTrigger>
      ),
    },
  ]);

  return (
    <Container fluid className="bg-white border rounded p-3">
      <LoadingOverlay active={isFetching}></LoadingOverlay>
      <h1 className="h4">Status</h1>
      <p>The status of your submitted refitting jobs will appear here.</p>
      <p>
        This list is tracked by your browser and will be lost if browser data is
        cleared.
      </p>
      <Button variant="light" className="mb-3" onClick={() => refetch()}>
        <i className="bi bi-arrow-clockwise" /> Refresh
      </Button>
      
      {/* Debug info */}
      <div className="mb-3">
        <small className="text-muted">
          Debug: {data.length} jobs found. 
          {data.length === 0 && " No jobs to display."}
        </small>
      </div>
      
      {data && data.length > 0 ? (
        <Table
          columns={columns}
          data={data}
          options={{
            initialState: {
              sortBy: [{ id: 'submittedAt', desc: true }],
            },
          }}
          striped
          bordered
        />
      ) : (
        <div className="text-center py-4">
          <p className="text-muted">No refitting jobs submitted yet.</p>
          <p className="text-muted">Submit a job to see it appear here.</p>
        </div>
      )}
    </Container>
  );
}
