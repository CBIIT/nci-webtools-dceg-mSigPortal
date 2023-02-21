import { useState, useMemo, useEffect } from 'react';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { Button } from 'react-bootstrap';
import { Link } from 'react-router-dom';
import { Container } from 'react-bootstrap';
import Table from '../../controls/table/table2';
import { useMultiJobStatusQuery } from './apiSlice';

export default function Status() {
  const [jobs, setJobs] = useState([]);
  const { data, isFetching, refetch } = useMultiJobStatusQuery(jobs);

  // get jobs from local storage
  useEffect(() => {
    const localJobs = JSON.parse(localStorage.getItem('jobs'));
    setJobs(localJobs);
  }, [localStorage.getItem('jobs')]);
  // update local storage
  useEffect(() => {
    if (jobs && jobs.length) localStorage.setItem('jobs', JSON.stringify(jobs));
  }, [jobs]);

  function removeJob(id) {
    const remaining = jobs.filter((e) => e !== id);
    localStorage.setItem('jobs', JSON.stringify(remaining));
    setJobs(remaining);
  }

  const columns = useMemo(() => [
    {
      accessor: 'jobName',
      Header: 'Name',
      Cell: (e) => (
        <Link exact to={`/extraction/${e.row.original.id}`}>
          {e.value}
        </Link>
      ),
    },
    {
      accessor: 'status',
      Header: 'Status',
      Cell: (e) => {
        if (e.value === 'SUBMITTED') return 'Submitted';
        else if (e.value === 'IN_PROGRESS') return 'In Progress';
        else if (e.value === 'FAILED') return 'Failed';
        else if (e.value === 'COMPLETED') return 'Completed';
        else return e.value;
      },
    },
    {
      Header: 'Remove',
      Cell: (e) => (
        <Button onClick={() => removeJob(e.row.original.id)}>Remove</Button>
      ),
      size: 100,
    },
  ]);

  return (
    <Container fluid className="bg-white border rounded p-3">
      <LoadingOverlay active={isFetching}></LoadingOverlay>
      <h4>Status</h4>
      <p>The status of your submitted jobs will appear here</p>
      <Button variant="light" className="mb-3" onClick={() => refetch()}>
        <i className="bi bi-arrow-clockwise" /> Refresh
      </Button>
      {data && jobs.length > 0 && (
        <Table className="border" columns={columns} data={data} striped />
      )}
    </Container>
  );
}
