import React from 'react';
import { Link } from 'react-router-dom';

export default function Instructions({ loading }) {
  const examples = [
    { title: 'VCF Example of User Input', path: 'vcfExample' },
    { title: 'Sherlock-Lung-232', path: 'sherlock-lung-232' },
    { title: 'Mutographs-ESCC', path: 'mutographs-escc' },
    { title: 'PCAWG Lung-AdenoCA', path: 'pcawg-LungAdenoCA' },
    { title: 'PCAWG Lung-SCC', path: 'pcawg-LungSCC' },
    { title: 'PCAWG Breast-AdenoCA', path: 'pcawg-BreastAdenoCA' },
    { title: 'PCAWG Skin-Melanoma', path: 'pcawg-SkinMelanoma' },
    { title: 'PCAWG PanCancer', path: 'pcawg-PanCancer' },
    { title: 'TCGA PanCancer', path: 'tcga-PanCancer' },
    // {
    //   title: 'MBD4 defect is associated with hypermutated CpG>TpG pattern',
    //   external: {
    //     name: 'PMID: 29760383',
    //     href: 'https://pubmed.ncbi.nlm.nih.gov/29760383/',
    //   },
    //   path: 'mbd4_defect',
    // },
  ];

  return (
    <div className="border rounded bg-white py-3 px-4">
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
