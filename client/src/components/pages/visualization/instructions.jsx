import { useSelector } from 'react-redux';
import { Table, Alert } from 'react-bootstrap';
import { useExampleHeaderQuery } from './userForm/apiSlice';
import { renderInlineCode } from '@/components/controls/utils/renderInlineCode';

export default function Instructions() {
  const { inputFormat } = useSelector((state) => state.visualization.userForm);
  const { source } = useSelector((state) => state.visualization.main);

  // get file header from example input files
  const { data, error } = useExampleHeaderQuery(
    inputFormat.value.toUpperCase(),
    {
      skip: !inputFormat?.value,
    }
  );

  const examples = [
    { title: 'VCF Example of User Input', path: 'vcfExample' },
    { title: 'Sherlock-Lung-232', path: 'sherlock-lung-232' },
    { title: 'Mutographs-ESCC', path: 'mutographs-escc' },
    { title: 'ChernobylThyroid', path: 'chernobyl-thyroid' },
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

  // structured column requirements per supported user input format
  const formatRequirements = {
    vcf: {
      description:
        'Upload a sample-level Variant Call Format file containing single-nucleotide variants and indels.',
      columns: [
        {
          column: '#CHROM',
          status: 'Required',
          requirement:
            'Chromosome containing the variant. Use one chromosome-naming style consistently.',
        },
        {
          column: 'POS',
          status: 'Required',
          requirement:
            'Positive, 1-based position in the selected reference genome.',
        },
        {
          column: 'ID',
          status: 'Required',
          requirement: 'Variant identifier; use "." when unavailable.',
        },
        { column: 'REF', status: 'Required', requirement: 'Reference allele.' },
        {
          column: 'ALT',
          status: 'Required',
          requirement:
            'One alternate allele. Multiallelic records with comma-separated ALT alleles are not supported.',
        },
        {
          column: 'QUAL',
          status: 'Required',
          requirement: 'Variant quality; use "." when unavailable.',
        },
        {
          column: 'FILTER',
          status: 'Required',
          requirement: 'Filter status or filter labels.',
        },
        {
          column: 'INFO',
          status: 'Required',
          requirement: 'VCF INFO field; use "." when unavailable.',
        },
        {
          column: 'FORMAT',
          status: 'Required',
          requirement:
            'Sample-field format. GT must be the first colon-delimited subfield.',
        },
        {
          column: 'Sample columns',
          status: 'Required; dynamic',
          requirement:
            'Include at least one sample column after FORMAT. Each column label is used as the sample identifier.',
        },
      ],
      fileFormat:
        'Tab-delimited VCF v4.x; use `.vcf` or `.vcf.gz`. `##` metadata lines may appear before the `#CHROM` header. Include one variant and one ALT allele per record.',
      example: '`demo_input_multi.vcf` or `demo_input_multi.vcf.gz`',
    },
    maf: {
      description:
        'Upload a sample-level Mutation Annotation Format table containing single-nucleotide variants and indels.',
      columns: [
        {
          column: 'Tumor_Sample_Barcode',
          status: 'Required',
          requirement: 'Non-empty sample identifier.',
        },
        {
          column: 'Chromosome',
          status: 'Required',
          requirement:
            'Chromosome containing the variant. Use one chromosome-naming style consistently.',
        },
        {
          column: 'Start_Position',
          status: 'Required',
          requirement:
            'Positive, 1-based start coordinate in the selected reference genome. The legacy spelling Start_position is accepted and mapped to Start_Position.',
        },
        {
          column: 'End_Position',
          status: 'Required',
          requirement:
            'Positive end coordinate using standard MAF conventions. The legacy spelling End_position is accepted and mapped to End_Position.',
        },
        {
          column: 'Reference_Allele',
          status: 'Required',
          requirement: 'Reference allele.',
        },
        {
          column: 'Tumor_Seq_Allele1',
          status: 'Required',
          requirement: 'First reported tumor allele.',
        },
        {
          column: 'Tumor_Seq_Allele2',
          status: 'Required',
          requirement:
            'Second reported tumor allele; the alternate allele is determined from the tumor and reference alleles.',
        },
      ],
      optionalColumns:
        'Additional standard MAF annotation fields are allowed and are not used to generate the mutational profile. Examples include `Hugo_Symbol`, `Entrez_Gene_Id`, `Center`, `NCBI_Build`, `Strand`, `Variant_Classification`, `Variant_Type`, `dbSNP_RS`, `dbSNP_Val_Status`, and `Matched_Norm_Sample_Barcode`.',
      fileFormat:
        'Tab-delimited text; use `.maf` or `.txt`. Include one variant per row. Coordinates and alleles should follow standard MAF representation and must match the reference genome selected in the form.',
      example: '`demo_input_multi_MAF.txt`',
    },
    csv: {
      description:
        'Upload a simple, row-based variant table. The file must contain exactly the seven columns below, in the order shown.',
      columns: [
        {
          column: 'SAMPLE',
          status: 'Required',
          requirement:
            'Sample identifier. Multiple rows may use the same sample identifier.',
        },
        {
          column: 'CHROM',
          status: 'Required',
          requirement: 'Chromosome containing the variant.',
        },
        {
          column: 'START',
          status: 'Required',
          requirement:
            'Positive start coordinate in the selected reference genome.',
        },
        {
          column: 'END',
          status: 'Required',
          requirement:
            'Positive end coordinate; it cannot be smaller than START.',
        },
        { column: 'REF', status: 'Required', requirement: 'Reference allele.' },
        {
          column: 'ALT',
          status: 'Required',
          requirement: 'One alternate allele.',
        },
        {
          column: 'FILTER',
          status: 'Required',
          requirement:
            'Filter status or label. Multiple filter labels in one cell may be separated by semicolons (;).',
        },
      ],
      fileFormat:
        'Comma-delimited text; use `.csv`. Do not add annotation columns, embedded commas, or blank lines. Include one variant and one ALT allele per row.',
      requiredHeader: '`SAMPLE,CHROM,START,END,REF,ALT,FILTER`',
      example: '`demo_input_multi.csv`',
    },
    tsv: {
      description:
        'Upload the same simple variant table used for CSV, with tabs as the delimiter. The file must contain exactly the seven columns below, in the order shown.',
      columns: [
        {
          column: 'SAMPLE',
          status: 'Required',
          requirement:
            'Sample identifier. Multiple rows may use the same sample identifier.',
        },
        {
          column: 'CHROM',
          status: 'Required',
          requirement: 'Chromosome containing the variant.',
        },
        {
          column: 'START',
          status: 'Required',
          requirement:
            'Positive start coordinate in the selected reference genome.',
        },
        {
          column: 'END',
          status: 'Required',
          requirement:
            'Positive end coordinate; it cannot be smaller than START.',
        },
        { column: 'REF', status: 'Required', requirement: 'Reference allele.' },
        {
          column: 'ALT',
          status: 'Required',
          requirement: 'One alternate allele.',
        },
        {
          column: 'FILTER',
          status: 'Required',
          requirement:
            'Filter status or label. Multiple filter labels in one cell may be separated by semicolons (;).',
        },
      ],
      fileFormat:
        'Tab-delimited text; use `.tsv`. Do not add annotation columns, embedded tabs, or blank lines. Include one variant and one ALT allele per row.',
      requiredHeader:
        '`SAMPLE<TAB>CHROM<TAB>START<TAB>END<TAB>REF<TAB>ALT<TAB>FILTER`',
      example: '`demo_input_multi.tsv`',
    },
    catalog_csv: {
      description:
        'Upload a precomputed mutational-count matrix. Rows are mutation categories and columns are samples.',
      columns: [
        {
          column: 'MutationType',
          status: 'Required',
          requirement:
            'First column. Each row must contain one unique mutation-category label from a single supported context.',
        },
        {
          column: 'Sample columns',
          status: 'Required; dynamic',
          requirement:
            'Include one or more uniquely named sample columns. Every cell must contain a finite, non-negative integer mutation count.',
        },
      ],
      fileFormat:
        'Comma-delimited text; use `.csv`. Do not include annotation or metadata columns, blank cells, `NA`, negative values, or nonnumeric values.',
      notes: [
        'The file must contain the complete mutation-category set for one supported context. Do not combine different contexts in one matrix. The detected context is shown during validation.',
      ],
      example: '`demo_input_catalog.csv`',
    },
    catalog_tsv: {
      description:
        'Upload a precomputed mutational-count matrix. Rows are mutation categories and columns are samples.',
      columns: [
        {
          column: 'MutationType',
          status: 'Required',
          requirement:
            'First column. Each row must contain one unique mutation-category label from a single supported context.',
        },
        {
          column: 'Sample columns',
          status: 'Required; dynamic',
          requirement:
            'Include one or more uniquely named sample columns. Every cell must contain a finite, non-negative integer mutation count.',
        },
      ],
      fileFormat:
        'Tab-delimited text; use `.tsv`. Do not include annotation or metadata columns, blank cells, `NA`, negative values, or nonnumeric values.',
      notes: [
        'The file must contain the complete mutation-category set for one supported context. Do not combine different contexts in one matrix. The detected context is shown during validation.',
      ],
      example: '`demo_input_catalog.tsv`',
    },
  };

  const requirements = formatRequirements[inputFormat?.value];

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
      {/* <hr />
      <h4>Example Queries</h4>
      <p>
        Choose an example query to view results for pre-selected parameters. You
        must reset between queries.
      </p>
      TBA */}
      {/* {examples.map(({ title, external, path }, index) => (
        <div key={index}>
          <Link to={`/visualization/example/${path}`}>
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
      {source == 'user' && (
        <>
          <hr />
          <h4>File requirements</h4>
          <p>
            Visualization accepts six primary formats: VCF, MAF, CSV, TSV,
            CATALOG CSV, and CATALOG TSV. VCF, MAF, CSV, and TSV contain
            individual variants. CATALOG CSV and CATALOG TSV contain a
            precomputed mutational-count matrix.
          </p>
          <p>
            For VCF, MAF, CSV, or TSV, select the reference genome build that
            matches the uploaded coordinates. Do not mix genome builds or
            chromosome-naming styles within a file. Reference-genome and
            experimental-strategy settings do not apply to catalog matrices.
          </p>
          <p>
            Under <b>Data Source: User</b>, choose a file format to view an
            example header. See the format-specific requirements below for
            required, optional, and sample-specific columns.
          </p>
          <Alert variant="info" className="mb-3">
            <ul
              className="mb-0 pl-3"
              style={{ display: 'block', listStyle: 'disc', columnCount: 1 }}
            >
              <li>
                Column names are case-insensitive. If a column name does not
                match the required name, an error is shown.
              </li>
              <li>
                Do not mix genome builds or chromosome-naming styles within a
                file.
              </li>
              <li>
                Reference-genome and experimental-strategy settings do not apply
                to catalog matrices.
              </li>
            </ul>
          </Alert>
          {requirements ? (
            <>
              <b>{inputFormat.label}</b>
              <p>{renderInlineCode(requirements.description)}</p>
              <Table striped bordered size="sm" responsive>
                <thead>
                  <tr>
                    <th>Column</th>
                    <th>Status</th>
                    <th>Requirement</th>
                  </tr>
                </thead>
                <tbody>
                  {requirements.columns.map((col) => (
                    <tr key={col.column}>
                      <td>{col.column}</td>
                      <td>{col.status}</td>
                      <td>{col.requirement}</td>
                    </tr>
                  ))}
                </tbody>
              </Table>
              {requirements.optionalColumns && (
                <p>
                  <b>Optional columns:</b>{' '}
                  {renderInlineCode(requirements.optionalColumns)}
                </p>
              )}
              <p>
                <b>File format:</b> {renderInlineCode(requirements.fileFormat)}
              </p>
              {requirements.requiredHeader && (
                <p>
                  <b>Required header:</b>{' '}
                  {renderInlineCode(requirements.requiredHeader)}
                </p>
              )}
              {requirements.notes?.map((note, index) => (
                <p key={index}>{renderInlineCode(note)}</p>
              ))}
              <p>
                <b>Example:</b> {renderInlineCode(requirements.example)} in{' '}
                <b>Download Example Data</b>.
              </p>
            </>
          ) : (
            <p>
              Select a file format under <b>Data Source: User</b> to view its
              column requirements.
            </p>
          )}
          <hr />
          <h4>Examples of file header for each supported format</h4>
          <p>
            Choose different file formats under <b>Data Source: User</b> to view
            different examples of file headers
          </p>
          <b>{inputFormat.label}</b>
          <pre className="border rounded bg-light p-3">{data}</pre>
        </>
      )}
    </div>
  );
}
