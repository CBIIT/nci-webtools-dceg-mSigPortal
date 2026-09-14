import { useEffect, useState } from 'react';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { Container, Table, Alert } from 'react-bootstrap';
import { useExampleQuery } from './apiSlice';
import { useHistory } from 'react-router-dom';
import { Button } from 'react-bootstrap';
import { renderInlineCode } from '../../controls/utils/renderInlineCode';

export default function Instructions({ formLimits }) {
  const examples = [
    {
      label: 'Example results based on SBS96 profiles',
      path: 'Example_SBS96_SigProfileExtractor',
    },
    {
      label: 'Example results based on DBS78 profiles',
      path: 'Example_DBS78_SigProfileExtractor',
    },
    {
      label: 'Example results based on ID83 profiles',
      path: 'Example_ID83_SigProfileExtractor',
    },
  ];
  const history = useHistory();
  const [id, setId] = useState(null);

  const { data: exampleData, isFetching } = useExampleQuery(id, {
    skip: !id,
  });

  useEffect(() => {
    if (exampleData && exampleData.id) {
      history.push(`/extraction/${exampleData.id}`);
    }
  }, [exampleData, history]);

  const handleExampleClick = (exampleFolder) => {
    setId(exampleFolder);
  };

  return (
    <Container fluid className="bg-white border rounded p-3">
      <LoadingOverlay active={isFetching} />
      <h4>Instructions</h4>
      <p>
        Choose a Data Source and its associated options to submit a query using
        the panel on the left
      </p>

      <hr />
      <div className="mt-2">
        <h4>Data Source</h4>
        <p>Public: Perform analysis using data available on the website</p>
        <p>User: Upload your own data</p>
        <div>
          <div>
            <a
              href="assets/exampleInput/extraction_sample_SBS96.all"
              download
              className="link-primary-underline"
            >
              Example of a SBS matrix
            </a>
          </div>
          <div>
            <a
              href="assets/exampleInput/extraction_sample_DBS78.all"
              download
              className="link-primary-underline"
            >
              Example of a DBS matrix
            </a>
          </div>
        </div>
      </div>
      <hr />
      <div className="mt-2 mb-2">
        <h4>Example Queries</h4>
        <p>
          Choose an example query to view results for pre-selected parameters.
          You must reset between queries.
        </p>
        <h6>SigProfilerExtraction</h6>

        {examples.map(({ label, path }, index) => (
          <Button
            key={index}
            onClick={() => handleExampleClick(path)}
            variant="link"
            className="d-block px-0"
          >
            {label}
          </Button>
        ))}
      </div>
      <hr />
      <div className="mt-2">
        <h4>User matrix file requirements</h4>
        <p>
          The Extraction user-data workflow accepts a mutational-count matrix,
          not raw VCF or MAF files. If starting from individual variants, first
          generate a matrix for the required mutation context and then upload
          that matrix.
        </p>
        <Alert variant="info" className="mb-3">
          Column names are case-insensitive. If a column name does not match the
          required name, an error is shown.
        </Alert>
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
                First column. Each row must contain one unique mutation-channel
                label that matches the selected Context Type. For example, SBS96
                requires the 96 standard trinucleotide channels and DBS78
                requires the 78 standard dinucleotide channels.
              </td>
            </tr>
            <tr>
              <td>Sample columns</td>
              <td>Required; dynamic</td>
              <td>
                Each remaining column represents one uniquely named sample.
                Values must be finite, non-negative integer mutation counts.
              </td>
            </tr>
          </tbody>
        </Table>
        <p>
          <b>File format:</b>
          {renderInlineCode(
            'Tab-delimited plain text. Recommended extensions are `.all`, `.tsv`, or `.txt`.'
          )}
        </p>
        <p className="mb-1">
          <b>Matrix rules:</b>
        </p>
        <ul style={{ display: 'block', listStyle: 'disc', columnCount: 1 }}>
          <li> Use one row per mutation channel and one column per sample.</li>
          <li>
            The mutation-channel names and number of rows must match the
            selected Context Type.
          </li>
          <li>
            {renderInlineCode(
              'Do not include duplicate mutation channels, duplicate sample names, blank cells, `NA`, negative values, or nonnumeric values.'
            )}
          </li>
          <li> Do not include an all-zero sample column.</li>
          <li>
            Select a reference genome, opportunity genome, context type, and
            reference-signature set that are compatible with the uploaded
            matrix.
          </li>
        </ul>
        <p>
          <b>Examples:</b>
          {renderInlineCode(
            '`extraction_sample_SBS96.all` for SBS96 and `extraction_sample_DBS78.all` for DBS78.'
          )}
        </p>
      </div>
      <hr />
      <div className="pt-2">
        <h6>SigProfilerExtraction Note</h6>
        <p>
          <b>Advanced Parameter Lower and Upper Limits</b>
        </p>
        <ol>
          <li>
            Minimum Signatures: ({formLimits.minimum_signatures[0]}-
            {formLimits.minimum_signatures[1]})
          </li>
          <li>
            Maximum Signatures: ({formLimits.maximum_signatures[0]}-
            {formLimits.maximum_signatures[1]})
          </li>
          <li>
            NMF Replicates Size: ({formLimits.nmf_replicates[0]}-
            {formLimits.nmf_replicates[1]})
          </li>
          <li>
            Minimum NMF Iterations: ({formLimits.min_nmf_iterations[0]}-
            {formLimits.min_nmf_iterations[1]})
          </li>
          <li>
            Maximum NMF Iterations: ({formLimits.max_nmf_iterations[0]}-
            {formLimits.max_nmf_iterations[1]})
          </li>
          <li>
            NMF Test Convergence: ({formLimits.nmf_test_conv[0]}-
            {formLimits.nmf_test_conv[1]})
          </li>
        </ol>
        <p>
          To enhance the performance of mutational signature extraction using
          SigProfilerExtractor, we recommend a two-step analysis approach with
          the following parameters:
        </p>
        <p>
          <b>Step 1: Determine the optimal number of detected signatures.</b>
        </p>

        <p>Update the Advanced Parameters (default values):</p>
        <ol>
          <li>Minimum Signatures: 1</li>
          <li>Maximum Signatures: 12</li>
          <li>NMF Replicates Size: 50</li>
          <li>Minimum NMF Iterations: 1000</li>
          <li>Maximum NMF Iterations: 100000</li>
          <li>NMF Test Convergence: 1000</li>
        </ol>
        <p>
          After completing this step, refer to the "Signature Map" table to
          identify the best number of detected signatures (denoted as N, for
          example).
        </p>
        <p>
          <b>Step 2: Refine the results of mutational signature extraction:</b>
        </p>
        <ol>
          <li>Minimum Signatures: N</li>
          <li>Maximum Signatures: N</li>
          <li>NMF Replicates Size: 100</li>
          <li>Minimum NMF Iterations: 10000</li>
          <li>Maximum NMF Iterations: 1000000</li>
          <li>NMF Test Convergence: 10000</li>
        </ol>
        <p>
          Performing these two-step analyses can significantly enhance the
          performance of mutational signature analysis using
          SigProfilerExtractor, particularly for jobs with large sample size.
        </p>
      </div>
    </Container>
  );
}
