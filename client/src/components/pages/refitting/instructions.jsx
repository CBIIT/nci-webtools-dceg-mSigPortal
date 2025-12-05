import React, { useEffect, useState } from 'react';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { Container } from 'react-bootstrap';
import { useHistory } from 'react-router-dom';
import { Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';

export default function Instructions() {
  const examples = [
    {
      label: 'Example results based on SBS96 targeted sequencing',
      path: 'Example_SBS96_Refitting',
    },
    {
      label: 'Example results based on DBS78 targeted sequencing',
      path: 'Example_DBS78_Refitting',
    },
  ];
  const history = useHistory();
  const [id, setId] = useState(null);
  const [isFetching, setIsFetching] = useState(false);
  
  // Get current signature type from Redux store
  const { signatureType } = useSelector((state) => state.refitting.userForm);

  const handleExampleClick = (exampleFolder) => {
    setIsFetching(true);
    setId(exampleFolder);
    // Simulate loading for example data
    setTimeout(() => {
      setIsFetching(false);
      // Navigate to refitting with example data loaded
    }, 1000);
  };

  return (
    <Container fluid className="bg-white border rounded p-3">
      <LoadingOverlay active={isFetching} />
      <h1 className='h4-title'>How to Submit a Query</h1>
      <p>
        Use the panel on the left to configure the required options.
      </p>

      <hr />
      
      <hr />
      <div className="mt-2">
        <h1 className='h4-title'>Signature Type</h1>
        <ul style={{ padding: 0, paddingLeft: '20px', columnCount: 'unset', columns: 'unset' }}>
          <li style={{ display: 'list-item', marginBottom: '8px' }}><strong>SBS:</strong> Single Base Substitutions</li>
          <li style={{ display: 'list-item', marginBottom: '8px' }}><strong>DBS:</strong> Double Base Substitutions</li>
        </ul>
      </div>
      <hr />
      <div className="mt-2">
        <h1 className='h4-title'>Reference Genome</h1>
        <ul style={{ padding: 0, paddingLeft: '20px', columnCount: 'unset', columns: 'unset' }}>
          <li style={{ display: 'list-item', marginBottom: '8px' }}><strong>hg19 (GRCh37):</strong> Genome Reference Consortium Human Build 37</li>
          <li style={{ display: 'list-item', marginBottom: '8px' }}><strong>hg38 (GRCh38):</strong> Genome Reference Consortium Human Build 38</li>
        </ul>
      </div>
      <hr />
      <div className="mt-2">
        <h1 className='h4-title'>Input Files</h1>
        <ul style={{ padding: 0, paddingLeft: '20px', columnCount: 'unset', columns: 'unset' }}>
          <li style={{ display: 'list-item', marginBottom: '8px' }}>
            <strong>MAF file:</strong> Contains {signatureType} mutation information for samples. <em>(Example provided)</em>
            <div>
            <a
              href={import.meta.env.BASE_URL + `/assets/examples/refitting/${signatureType}_MAF_two_samples.txt`}
              download
              className="link-primary-underline"
            >
              Example of a {signatureType} MAF file
            </a>
          </div>
          </li>
          <li style={{ display: 'list-item', marginBottom: '8px' }}>
            <strong>Genomic file:</strong> Defines the genomic regions targeted by sequencing panels. <em>(Example provided)</em>
          <div>
            <a
              href={import.meta.env.BASE_URL + "/assets/examples/refitting/Genomic_information_sample.txt"}
              download
              className="link-primary-underline"
            >
              Example of a genomic file
            </a>
          </div>
          </li>
          <li style={{ display: 'list-item', marginBottom: '8px' }}>
            <strong>Clinical file:</strong> Specifies sample ID, sequencing panel ID, and cancer type.
            <ul style={{ padding: 0, paddingLeft: '20px', marginTop: '4px', columnCount: 'unset', columns: 'unset' }}>
              <li style={{ display: 'list-item' }}>Cancer type must match one from the cancer type dictionary file.</li>
            </ul>
            <div>
            <a
              href={import.meta.env.BASE_URL + `/assets/examples/refitting/${signatureType}Clinical_sample.txt`}
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
        <h1 className='h4-title'>Output File</h1>
        <ul style={{ padding: 0, paddingLeft: '20px', columnCount: 'unset', columns: 'unset' }}>
          <li style={{ display: 'list-item', marginBottom: '8px' }}>
            <strong>Signature activity and burden:</strong> Provides both the estimated activity of each signature 
            and the estimated number of mutations caused by each signature for each subject.
          </li>
        </ul>
      </div>
      
    </Container>
  );
}
