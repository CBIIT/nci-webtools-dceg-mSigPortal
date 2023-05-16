import React from 'react';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { Link } from 'react-router-dom';
import { Container } from 'react-bootstrap';
import { useExampleQuery, useRefreshQuery } from './apiSlice';
import { useLocation } from 'react-router-dom';
import { v4 as uuidv4 } from 'uuid';

export default function Instructions({ props, loading }) {
  const examples = [
    {
      title: 'Example results based on SBS96 profiles',
      path: 'Example_SBS96_SigProfileExtractor',
    },
    {
      title: 'Example results based on DBS78 profiles',
      path: 'Example_DBS78_SigProfileExtractor',
    },
    {
      title: 'Example results based on ID83 profiles',
      path: 'Example_ID83_SigProfileExtractor',
    },
  ];
  const location = useLocation();
  const pathParts = location.pathname.split('/');
  const id = pathParts[pathParts.length - 1];
  const randomUUID = uuidv4();

  function ExampleComponent({ id }) {
    const { data: exampleData, refetch: refresh } = useExampleQuery(id, {
      skip: !id,
    });
  }
  if (id.startsWith('Example_')) {
    ExampleComponent({ id });
  }

  // const { data: exampleData, refetch: refresh } = useExampleQuery(id, {
  //   skip: !id,
  // });

  return (
    <Container fluid className="bg-white border rounded p-3" {...props}>
      <LoadingOverlay active={loading} />
      <h4>Instructions</h4>
      <p>
        Choose a Data Source and its associated options to submit a query using
        the panel on the left
      </p>
      <hr />
      <h4>Data Source</h4>
      <p>Public: Perform analysis using data available on the website</p>
      <p>User: Upload your own data</p>
      <hr />
      <h4>Example Queries</h4>
      <p>
        Choose an example query to view results for pre-selected parameters. You
        must reset between queries.
      </p>
      <h6>SigProfilerExtraction</h6>

      {examples.map(({ title, external, path }, index) => (
        <div key={index}>
          <Link to={`/extraction/${path}_${randomUUID}`} disabled>
            <span className="sr-only">{title + ' link'}</span>
            {title}
          </Link>
          {external && (
            <span>
              {'; '}
              <a href={external.href} target="_blank">
                {external.name}
              </a>
            </span>
          )}
        </div>
      ))}
    </Container>
  );
}
