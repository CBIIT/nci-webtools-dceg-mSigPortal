import React from 'react';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { Link } from 'react-router-dom';
import { Container } from 'react-bootstrap';

export default function Instructions({ props, loading }) {
  const examples = [
    // {
    //   title: 'Sherlock-Lung-232',
    //   path: 'sherlock-lung-232',
    // },
    // {
    //   title: 'Mutographs-ESCC',
    //   path: 'mutographs-escc',
    // },
    // {
    //   title: 'PCAWG Lung-AdenoCA',
    //   path: 'pcawg-lungadenoca',
    // },
    // {
    //   title: 'PCAWG Lung-SCC',
    //   path: 'pcawg-lungscc',
    // },
    // {
    //   title: 'PCAWG Breast-AdenoCA',
    //   path: 'pcawg-breastadenoca',
    // },
    // {
    //   title: 'PCAWG Skin-Melanoma',
    //   path: 'pcawg-skinmelanoma',
    // },
  ];

  const exampleFiles = [
    {
      title: 'Example results based on SBS96 profiles',
      path: 'extraction-SBS96',
    },
    {
      title: 'Example results based on DBS78 profiles',
      path: 'extraction-DBS78',
    },
    {
      title: 'Example results based on ID83 profiles',
      path: 'extraction-ID83',
    },
  ];

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
      <h5>SigProfilerExtraction</h5>

      {exampleFiles.map(({ title, external, path }, index) => (
        <div key={index}>
          <Link to={`/extraction/${path}`} disabled>
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
