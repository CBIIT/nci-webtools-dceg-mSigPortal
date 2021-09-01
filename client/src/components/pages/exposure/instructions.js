import React from 'react';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { Link } from 'react-router-dom';

export default function Instructions({ loading }) {
  const examples = [
    {
      name: 'Lung',
      title:
        'PCAWG/WGS/COSMIC_v3_Signatures_GRCh37_SBS96/ Lung-AdenoCA; MSA SBS5 vs SBS40',
      path: 'pcawg-lungadenoca',
    },
    {
      name: 'Skin',
      title:
        'PCAWG/WGS/COSMIC_v3_Signatures_GRCh37_SBS96/ Skin-Melanoma; MSA SBS7a vs SBS7b',
      path: 'pcawg-skinmelanoma',
    },
    {
      name: 'Breast',
      title:
        'PCAWG/WGS/COSMIC_v3_Signatures_GRCh37_SBS96/ Breast-AdenoCA; MSA SBS3 vs SBS5',
      path: 'pcawg-breastadenoca',
    },
  ];

  return (
    <div className="py-3 px-4">
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

      {examples.map(({ title, external, path }, index) => (
        <div key={index}>
          <Link to={`/exploration/${path}`} disabled>
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
    </div>
  );
}
