import { useState, useMemo, useEffect } from 'react';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { Button, Popover, OverlayTrigger } from 'react-bootstrap';
import { Link } from 'react-router-dom';
import { Container } from 'react-bootstrap';
import Table from '../../controls/table/table2';
import { useMultiJobStatusQuery } from './apiSlice';
import moment from 'moment';
import momentDurationFormatSetup from 'moment-duration-format';
momentDurationFormatSetup(moment);

export default function Status() {
  const [jobs, setJobs] = useState([]);
  const { data, isFetching, refetch } = useMultiJobStatusQuery(jobs, {
    skip: !jobs || !jobs.length,
  });

  // get jobs from local storage
  useEffect(() => {
    const localJobs = JSON.parse(localStorage.getItem('refitting-jobs') || '[]');
    setJobs(localJobs);
  }, []);

  // update jobs in local storage
  useEffect(() => {
    if (jobs && jobs.length) localStorage.setItem('refitting-jobs', JSON.stringify(jobs));
  }, [jobs]);

  function removeJob(id) {
    const remaining = jobs.filter((e) => e !== id);
    setJobs(remaining);
  }

  const columns = useMemo(() => [
    {
      accessor: 'jobName',
      Header: 'Name',
      Cell: ({ row, value }) => (
        <Link exact to={`/refitting/${row.original.id}`}>
          {value}
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
        else return value;
      },
    },
    {
      accessor: 'submittedAt',
      Header: 'Submitted Time',
      Cell: ({ value }) => {
        const time = moment(value);
        return `${time.format('LLL')} (${time.fromNow()})`;
      },
    },
    {
      id: 'duration',
      accessor: 'stopped',
      Header: 'Duration',
      Cell: ({ row, value }) => {
        if (value) {
          const start = moment(row.original.submittedAt);
          const end = moment(value);
          return moment
            .duration(end.diff(start), 'milliseconds')
            .format('hh:mm:ss', { trim: false });
        } else return '';
      },
    },
    {
      id: 'Download',
      accessor: 'status',
      Header: 'Download',
      Cell: ({ row, value }) => (
        <Button
          variant="link"
          href={`api/downloadRefittingOutput/${row.original.id}`}
          disabled={!['COMPLETED', 'FAILED'].includes(value)}
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
      {data && data.length > 0 && (
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
      )}
    </Container>
  );
}
