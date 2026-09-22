import { LoadingOverlay } from '@/components/controls/loading-overlay/loading-overlay';
import { Container, Table, Alert } from 'react-bootstrap';
import { renderInlineCode } from '@/components/controls/utils/renderInlineCode';

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
      <div className="mt-2">
        <h4>User data file requirements</h4>
        <p>
          Exploration requires an exposure/activity file and a mutational-count
          matrix. A signature-profile file is also required when{' '}
          <b>Use Public Signature Data</b> is turned off. All primary files must
          describe the same mutation context.
        </p>
        <Alert variant="info" className="mb-3">
          Use tab-delimited plain-text files with{' '}
          {renderInlineCode('`.tsv` or `.txt`')} extensions. Column names are
          case-insensitive; if a column name does not match the required name,
          an error is shown.
        </Alert>

        <h6>Exposure/activity file</h6>
        <Table striped bordered size="sm" responsive>
          <thead>
            <tr>
              <th>Column</th>
              <th>Status</th>
              <th>Requirement</th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td>Samples</td>
              <td>Required</td>
              <td>
                First column. Include one row per unique sample identifier.
              </td>
            </tr>
            <tr>
              <td>Signature columns</td>
              <td>Required; dynamic</td>
              <td>
                Each remaining column is one unique signature name. Values must
                be finite, non-negative integer mutation activities.
              </td>
            </tr>
          </tbody>
        </Table>
        <p>
          <b>Example:</b> {renderInlineCode('`Sherlock_SBS96_exposure.txt`')}.
        </p>

        <h6>Mutational-count matrix</h6>
        <Table striped bordered size="sm" responsive>
          <thead>
            <tr>
              <th>Column</th>
              <th>Status</th>
              <th>Requirement</th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td>MutationType</td>
              <td>Required</td>
              <td>
                First column. Include one unique mutation channel per row.
              </td>
            </tr>
            <tr>
              <td>Sample columns</td>
              <td>Required; dynamic</td>
              <td>
                Each remaining column is one unique sample identifier. Values
                must be finite, non-negative integer mutation counts.
              </td>
            </tr>
          </tbody>
        </Table>
        <p>
          <b>Example:</b> {renderInlineCode('`Sherlock_SBS96_matrix.txt`')}.
        </p>

        <h6>Signature-profile file</h6>
        <p className="text-muted mb-2">
          Required only when <b>Use Public Signature Data</b> is turned off.
        </p>
        <Table striped bordered size="sm" responsive>
          <thead>
            <tr>
              <th>Column</th>
              <th>Status</th>
              <th>Requirement</th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td>MutationType</td>
              <td>Required</td>
              <td>
                First column. Mutation-channel labels must match the
                mutational-count matrix exactly.
              </td>
            </tr>
            <tr>
              <td>Signature columns</td>
              <td>Required; dynamic</td>
              <td>
                Each remaining column is one unique signature name. Values must
                be finite and non-negative. Each signature column should sum to
                1 within numerical tolerance.
              </td>
            </tr>
          </tbody>
        </Table>
        <p>
          <b>Example:</b> {renderInlineCode('`Sherlock_SBS96_signature.txt`')}.
        </p>

        <p className="mb-1">
          <b>Cross-file compatibility:</b>
        </p>
        <ul style={{ display: 'block', listStyle: 'disc', columnCount: 1 }}>
          <li>
            Sample IDs in the exposure file should match the matrix
            sample-column labels.
          </li>
          <li>
            Signature names in the exposure file should match the selected or
            uploaded signature profiles.
          </li>
          <li>
            The matrix and signature-profile files should contain the same
            mutation-channel set, and all files should represent the same
            context (such as SBS96).
          </li>
          <li>
            Select the genome and public signature set that correspond to that
            context.
          </li>
          <li>
            Samples, signatures, or mutation channels that do not match across
            files are ignored rather than reconciled, so align them before
            uploading.
          </li>
        </ul>
      </div>
      {/* <hr />
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
      ))} */}
    </Container>
  );
}
