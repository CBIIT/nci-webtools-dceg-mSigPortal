import { Container, Table } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { renderInlineCode } from '../../controls/utils/renderInlineCode';

export default function Instructions() {
  const { signatureType } = useSelector((state) => state.refitting.userForm);

  return (
    <Container fluid className="bg-white border rounded p-3">
      <h1 className="h4-title">Refitting for Targeted Sequencing</h1>
      <p>
        This webpage provides mutational signature refitting for targeted
        sequencing data. The set of mutational signatures used for refitting for
        each cancer type is shown on the{' '}
        <a
          href={`#/catalog/sts`}
          className="link-primary-underline"
          target="_blank"
        >
          signature catalogue webpage
        </a>
        . The applicable cancer types must be specified in the clinical input
        file and have to be matched to those listed in the Cancer Type column of
        the{' '}
        <a
          href="assets/examples/refitting/CancerTypes_Dictionary.csv"
          download
          className="link-primary-underline"
        >
          cancer dictionary file
        </a>
        . Use the panel on the left to configure the remaining required inputs.
      </p>
      <hr />
      <div className="mt-2">
        <h1 className="h4-title">Signature Type</h1>
        <ul
          style={{
            padding: 0,
            paddingLeft: '20px',
            columnCount: 'unset',
            columns: 'unset',
          }}
        >
          <li style={{ display: 'list-item', marginBottom: '8px' }}>
            <strong>SBS:</strong> Single Base Substitutions
          </li>
          <li style={{ display: 'list-item', marginBottom: '8px' }}>
            <strong>DBS:</strong> Double Base Substitutions
          </li>
        </ul>
      </div>
      <hr />
      <div className="mt-2">
        <h1 className="h4-title">Reference Genome</h1>
        <ul
          style={{
            padding: 0,
            paddingLeft: '20px',
            columnCount: 'unset',
            columns: 'unset',
          }}
        >
          <li style={{ display: 'list-item', marginBottom: '8px' }}>
            <strong>hg19 (GRCh37):</strong> Genome Reference Consortium Human
            Build 37
          </li>
          <li style={{ display: 'list-item', marginBottom: '8px' }}>
            <strong>hg38 (GRCh38):</strong> Genome Reference Consortium Human
            Build 38
          </li>
        </ul>
      </div>
      <hr />
      <div className="mt-2">
        <h1 className="h4-title">Input Files</h1>
        <ul
          style={{
            padding: 0,
            paddingLeft: '20px',
            columnCount: 'unset',
            columns: 'unset',
          }}
        >
          <li style={{ display: 'list-item', marginBottom: '8px' }}>
            <strong>MAF file:</strong> Contains {signatureType} mutation
            information for samples. <em>(Example provided)</em>
            <div>
              <a
                href={`assets/examples/refitting/${signatureType}_MAF_two_samples.txt`}
                download
                className="link-primary-underline"
              >
                Example of a {signatureType} MAF file
              </a>
            </div>
          </li>
          <li style={{ display: 'list-item', marginBottom: '8px' }}>
            <strong>Genomic file:</strong> Defines the genomic regions targeted
            by sequencing panels. <em>(Example provided)</em>
            <div>
              <a
                href={
                  'assets/examples/refitting/Genomic_information_sample.txt'
                }
                download
                className="link-primary-underline"
              >
                Example of a genomic file
              </a>
            </div>
          </li>
          <li style={{ display: 'list-item', marginBottom: '8px' }}>
            <strong>Clinical file:</strong> Specifies sample ID, sequencing
            panel ID, and cancer type.
            <ul
              style={{
                padding: 0,
                paddingLeft: '20px',
                marginTop: '4px',
                columnCount: 'unset',
                columns: 'unset',
              }}
            >
              <li>
                {' '}
                Cancer type must match one from the cancer type dictionary file.
              </li>
            </ul>
            <div>
              <a
                href={`assets/examples/refitting/${signatureType}_Clinical_sample.txt`}
                download
                className="link-primary-underline"
              >
                Clinical sample file
              </a>
            </div>
          </li>
        </ul>
      </div>
      <hr />
      <div className="mt-2">
        <h1 className="h4-title">Input file columns</h1>
        <h6>Mutation MAF file</h6>
        <Table striped bordered size="sm" responsive>
          <thead>
            <tr>
              <th>Column</th>
              <th>Requirement</th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td>Chromosome</td>
              <td>Use 1–22, X, or Y. Do not add a chr prefix.</td>
            </tr>
            <tr>
              <td>Start_Position</td>
              <td>
                Positive, 1-based coordinate. Capitalization variants such as
                Start_position are mapped to Start_Position.
              </td>
            </tr>
            <tr>
              <td>End_Position</td>
              <td>
                For SBS, use the same coordinate as Start_Position. For DBS, use
                the coordinate of the second adjacent base.
              </td>
            </tr>
            <tr>
              <td>Variant_Type</td>
              <td>
                Use SNP for SBS and DNP for DBS. Values are case-sensitive. Do
                not mix SBS and DBS rows in one file.
              </td>
            </tr>
            <tr>
              <td>Reference_Allele</td>
              <td>
                Reference allele: one A/C/G/T base for SBS or two adjacent
                A/C/G/T bases for DBS.
              </td>
            </tr>
            <tr>
              <td>Tumor_Seq_Allele2</td>
              <td>
                Alternate allele: one A/C/G/T base for SBS or two adjacent
                A/C/G/T bases for DBS.
              </td>
            </tr>
            <tr>
              <td>Tumor_Sample_Barcode</td>
              <td>
                Non-empty sample identifier. It must match one SAMPLE_ID in the
                clinical file.
              </td>
            </tr>
          </tbody>
        </Table>
        <p>
          Additional MAF annotation columns are allowed and are ignored by
          Refitting. <b>File format:</b>{' '}
          {renderInlineCode(
            'tab-delimited `.txt` or `.maf`, one mutation per row.'
          )}
        </p>

        <h6>Genomic/panel-region file</h6>
        <Table striped bordered size="sm" responsive>
          <thead>
            <tr>
              <th>Column</th>
              <th>Requirement</th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td>Chromosome</td>
              <td>Use 1–22, X, or Y. Do not add a chr prefix.</td>
            </tr>
            <tr>
              <td>Start_Position</td>
              <td>
                Positive, 1-based start coordinate of the targeted interval.
              </td>
            </tr>
            <tr>
              <td>End_Position</td>
              <td>
                Inclusive end coordinate. It must be equal to or greater than
                Start_Position.
              </td>
            </tr>
            <tr>
              <td>SEQ_ASSAY_ID</td>
              <td>
                Sequencing-panel identifier. Multiple intervals may use the same
                panel identifier.
              </td>
            </tr>
          </tbody>
        </Table>
        <p>
          <b>File format:</b>{' '}
          {renderInlineCode('tab-delimited `.txt`, one interval per row.')}
        </p>

        <h6>Clinical file</h6>
        <Table striped bordered size="sm" responsive>
          <thead>
            <tr>
              <th>Column</th>
              <th>Requirement</th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td>SAMPLE_ID</td>
              <td>
                Unique sample identifier. It must match Tumor_Sample_Barcode in
                the mutation MAF. Use one clinical row per sample.
              </td>
            </tr>
            <tr>
              <td>SEQ_ASSAY_ID</td>
              <td>
                Sequencing-panel identifier. It must match a panel identifier in
                the genomic/panel-region file.
              </td>
            </tr>
            <tr>
              <td>CANCER_TYPE</td>
              <td>
                Cancer type. It must match a CANCER_TYPE in the downloadable
                cancer-type dictionary.
              </td>
            </tr>
          </tbody>
        </Table>
        <p>
          <b>File format:</b>{' '}
          {renderInlineCode('tab-delimited `.txt`, one row per sample.')}
        </p>

        <p className="mb-1">
          <b>Required matching across files:</b>
        </p>
        <ul style={{ display: 'block', listStyle: 'disc', columnCount: 1 }}>
          <li>
            {renderInlineCode(
              'MAF `Tumor_Sample_Barcode` ↔ clinical `SAMPLE_ID`.'
            )}
          </li>
          <li>
            {renderInlineCode(
              'Clinical `SEQ_ASSAY_ID` ↔ genomic-file `SEQ_ASSAY_ID`.'
            )}
          </li>
          <li>
            {renderInlineCode(
              'Clinical `CANCER_TYPE` ↔ cancer-dictionary `CANCER_TYPE`.'
            )}
          </li>
          <li>
            These identifiers are matched by value and are case-sensitive.
            Samples without a matching clinical record, or panels without a
            genomic interval, are skipped rather than analyzed, so align them
            before submission.
          </li>
        </ul>
      </div>
      <hr />
      <div className="mt-2">
        <h1 className="h4-title">Output File</h1>
        <ul
          style={{
            padding: 0,
            paddingLeft: '20px',
            columnCount: 'unset',
            columns: 'unset',
          }}
        >
          <li style={{ display: 'list-item', marginBottom: '8px' }}>
            <strong>Signature activity and burden:</strong> Provides both the
            estimated activity of each signature and the estimated number of
            mutations caused by each signature for each subject.
          </li>
        </ul>
      </div>
    </Container>
  );
}
