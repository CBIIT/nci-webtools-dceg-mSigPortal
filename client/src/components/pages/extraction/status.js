import { useState, useMemo, useEffect } from 'react';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { Button, Popover, OverlayTrigger } from 'react-bootstrap';
import { Link } from 'react-router-dom';
import { Container } from 'react-bootstrap';
import Table from '../../controls/table/table2';
import { useMultiJobStatusQuery } from './apiSlice';
import moment from 'moment';

export default function Status() {
  const [jobs, setJobs] = useState([]);
  const { data, isFetching, refetch } = useMultiJobStatusQuery(jobs, {
    skip: !jobs || !jobs.length,
  });

  // get jobs from local storage
  useEffect(() => {
    const localJobs = JSON.parse(localStorage.getItem('jobs'));
    setJobs(localJobs);
  }, [localStorage.getItem('jobs')]);
  // update jobs in local storage
  useEffect(() => {
    if (jobs && jobs.length) localStorage.setItem('jobs', JSON.stringify(jobs));
  }, [jobs]);

  function removeJob(id) {
    const remaining = jobs.filter((e) => e !== id);
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
      accessor: 'submittedAt',
      Header: 'Submitted Time',
      Cell: (e) => {
        const time = moment(e.value);
        return `${time.format('LLL')} (${time.fromNow()})`;
      },
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
      <h4>Status</h4>
      <p>The status of your submitted jobs will appear here.</p>
      <p>
        This list is tracked by your browser and will be lost if browser data is
        cleared.
      </p>
      <Button variant="light" className="mb-3" onClick={() => refetch()}>
        <i className="bi bi-arrow-clockwise" /> Refresh
      </Button>
      {data && data.length > 0 && (
        <Table columns={columns} data={data} striped bordered />
      )}
    </Container>
  );
}
