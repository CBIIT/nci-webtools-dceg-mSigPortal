import React from 'react';
import { Link } from 'react-router-dom';

export default function Instructions({ loading }) {
  const examples = [
    { title: 'VCF Example', path: 'vcfExample' },
    { title: 'PCAWG Lung-AdenoCA', path: 'pcawg-lungadenoca' },
    { title: 'PCAWG Lung-SCC', path: 'pcawg-lungscc' },
    { title: 'PCAWG Breast-AdenoCA', path: 'pcawg-breastadenoca' },
    { title: 'PCAWG Skin-Melanoma', path: 'pcawg-skinmelanoma' },
    { title: 'PCAWG PanCancer', path: 'pcawg-pancancer' },
    { title: 'TCGA PanCancer', path: 'tcga-pancancer' },
    {
      title: 'MBD4 defect is associated with hypermutated CpG>TpG pattern',
      external: {
        name: 'PMID: 29760383',
        href: 'https://pubmed.ncbi.nlm.nih.gov/29760383/',
      },
      path: 'mbd4_defect',
    },
  ];

  return (
    <div className="border rounded bg-white p-3">
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
      {examples.map(({ title, external, path }, index) => (
        <div key={index}>
          <Link to={`/visualization/example/${path}`} disabled>
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
