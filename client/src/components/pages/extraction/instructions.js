import React, { useEffect, useState, useCallback } from 'react';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { Link } from 'react-router-dom';
import { Container } from 'react-bootstrap';
import { useExampleQuery, useRefreshQuery } from './apiSlice';
import { useHistory } from 'react-router-dom';
import { Button } from 'react-bootstrap';

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
  const history = useHistory();

  //const id = pathParts[pathParts.length - 1];

  const [uuid, setUuid] = useState(null);

  const [id, setId] = useState(null);

  const { data: exampleData } = useExampleQuery(id, {
    skip: !id,
  });

  console.log('exampleData', exampleData);

  useEffect(() => {
    if (exampleData && exampleData.id) {
      history.push(`/extraction/${exampleData.id}`);
    }
  }, [exampleData, history]);

  const handleExampleClick = async (exampleFolder) => {
    try {
      // Fetch example data here if needed
      setId(exampleFolder);
    } catch (error) {
      console.error(error); // Handle the error
    }
  };

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
        <div>
          <Button
            key={index}
            onClick={() => handleExampleClick(path)}
            variant="link"
          >
            <span className="sr-only">{title + ' link'}</span>
            {title}
            {external && (
              <span>
                {'; '}
                <a href={external.href} target="_blank">
                  {external.name}
                </a>
              </span>
            )}
          </Button>
        </div>
      ))}
    </Container>
  );
}
