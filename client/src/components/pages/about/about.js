import React from 'react';
import { Container } from 'react-bootstrap';

export default function About() {
  return (
    <div className="mx-3">
      <div className="bg-white border p-3 mx-3" style={{ minHeight: '500px' }}>
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

          <div id="catalog">
            <h5 class="text-catalog">Signature Catalog</h5>
            <p>
              Comprehensively exploring curated census of mutational signatures
              from scientific literature, currently including COSMIC mutational
              signatures, environmental mutagenesis, DNA repair gene edits,
              cancer specific signatures and others. In addition, it allows
              users to visualize, compare and download these reference
              signatures.
            </p>
          </div>
          <div id="visualization">
            <h5 class="text-visualization">Signature Visualization</h5>
            <p>
              Interactively visualizing and analyzing mutational profiles at the
              sample level from both user input (formats supported include VCF,
              MAF, CSV, TXT, catalog, etc) and collected cancer genomic studies
              from scientific literature. It allows users to perform a wide
              range of analyses including cosine similarity, enrichment
              analysis, mutational profile comparison, and principal components
              analysis for all different types of mutational profiles including
              SBS, INDEL, DBS, SV, and CNV. Kataegis identification is also
              supported for VCF input format. With respect to user input, all
              visualizations and analyses can be performed across different
              provided groups of assigned mutations in the same sample (such as
              clustered mutations vs. non-clustered mutations; nonsynonymous
              mutations vs. other mutations, etc).
            </p>
          </div>
          <div id="extraction">
            <h5 class="text-extraction">Signature Extraction</h5>
            <p>
              This module allows users to perform both mutational signature de
              novo extraction and decomposition analysis using state-of-the-art
              algorithms such as SigProfiler, Signal, MuSiCal, and
              SignatureAnalyzer. It supports both user input and mutation
              profile data collected from popular cancer genomic studies. Users
              can select any combination of known mutational signatures for the
              signature decomposition analysis as a reference. The results from
              this module can be fully imported into the Signature Exploration
              module for visualization and comparison.
            </p>
          </div>
          <div id="exploration">
            <h5 class="text-exploration">Signature Exploration</h5>
            <p>
              Systematically exploring the mutational signature activities and
              performance of the mutational signature decomposition from user
              input, or collected public genomic studies (
              <a
                href="https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga"
                target="_blank"
                rel="noreferrer"
              >
                TCGA
              </a>
              ,{' '}
              <a
                href="https://dcc.icgc.org/pcawg"
                target="_blank"
                rel="noreferrer"
              >
                PCAWG
              </a>
              ,{' '}
              <a
                href="https://dceg.cancer.gov/research/cancer-types/lung/sherlock-lung-study"
                target="_blank"
                rel="noreferrer"
              >
                Sherlock-Lung
              </a>
              , etc.). This module allows users to perform analyses with
              mutational signature patterns and integratively explore the
              activities of each mutational signature including visualizations
              of tumor mutational burden, signature decomposition performance,
              mutational signature associations, sample clustering by mutational
              signatures, prevalence of single mutational signatures, and
              decomposition of mutational signatures in individual samples.
            </p>
          </div>
          <div id="association">
            <h5 class="text-association">Signature Association</h5>
            <p>
              Statistically analyzing and visualizing associations between
              mutational signature activities (using different measurements) and
              collected sample level variables including genomic features,
              epigenomic features, mutational status, copy number alterations or
              clinical variables from different cancer genomic studies. In
              addition, this module allows users to select different statistical
              approaches for both univariable and multivariable association
              analyses.
            </p>
          </div>
          <div id="api">
            <h5 class="text-api">Signature API Access</h5>
            <p>
              The mutational signature-related data in mSigPortal can be
              accessed programmatically through a REST API. The syntax for
              accessing the data is similar to the web address link used for
              queries on the webpage for each mSigPortal module, and the text
              output returned is typically the same as the file that can be
              downloaded from the online site. Users can connect directly to the
              API to perform batch queries and develop their own visualizations
              or analyses.
            </p>
          </div>
          <div>
            {' '}
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
            <br />
            <div>
              <h5>Release History</h5>
              <i>mSigPortal 1.0.0</i>
              <ul>
                <li>Initial Release</li>
              </ul>
            </div>
          </div>
        </Container>
      </div>
    </div>
  );
}
