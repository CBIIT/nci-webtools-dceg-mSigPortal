import React from 'react';
import { Container } from 'react-bootstrap';
import { useSelector } from 'react-redux';

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
          href="/#/catalog/sts"
          className="link-primary-underline"
          target="_blank"
        >
          signature catalogue webpage
        </a>
        . The applicable cancer types must be specified in the clinical input
        file and have to be matched to those listed in the Cancer Type column of
        the{' '}
        <a
          href="/assets/examples/refitting/CancerTypes_Dictionary.csv"
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
                href={`/assets/examples/refitting/${signatureType}_MAF_two_samples.txt`}
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
                  '/assets/examples/refitting/Genomic_information_sample.txt'
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
              <li style={{ display: 'list-item' }}>
                Cancer type must match one from the cancer type dictionary file.
              </li>
            </ul>
            <div>
              <a
                href={`/assets/examples/refitting/${signatureType}_Clinical_sample.txt`}
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
