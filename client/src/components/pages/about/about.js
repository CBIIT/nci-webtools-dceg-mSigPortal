import React from 'react';
import { Container } from 'react-bootstrap';

export default function About() {
  return (
    <div className="mx-3">
      <div className="bg-white border p-3 mx-3" style={{ minHeight: '420px' }}>
        <Container fluid="xl">
          <p>
            mSigPortal was developed by Dr.{' '}
            <a
              href="https://dceg.cancer.gov/about/staff-directory/zhang-tongwu"
              target="_blank"
              rel="noreferrer"
            >
              Tongwu Zhang
            </a>{' '}
            in the lab of Dr.{' '}
            <a
              href="https://dceg.cancer.gov/about/staff-directory/landi-maria"
              target="_blank"
              rel="noreferrer"
            >
              Maria Teresa Landi
            </a>{' '}
            in collaboration with the NCI Center for Biomedical Informatics and
            Information Technology (CBIIT). Support for this project comes from
            the Division of Cancer Epidemiology and Genetics Informatics Tool
            Challenge as well as the{' '}
            <a
              href="https://dceg.cancer.gov/about/organization/tdrp/iteb"
              target="_blank"
              rel="noreferrer"
            >
              Integrative Tumor Epidemiology Branch
            </a>
            .
          </p>
          <p>
            mSigPortal’s{' '}
            <a
              href="https://github.com/CBIIT/nci-webtools-dceg-mSigPortal"
              target="_blank"
              rel="noreferrer"
            >
              source code
            </a>{' '}
            is available in Github under MIT license, an Open Source Initiative
            approved license.
          </p>
          <p>
            Questions or comments? Contact{' '}
            <a
              href="mailto:NCImSigPortalWebAdmin@mail.nih.gov"
              target="_blank"
              rel="noreferrer"
            >
              support
            </a>
            .
          </p>
        </Container>
      </div>
    </div>
  );
}
