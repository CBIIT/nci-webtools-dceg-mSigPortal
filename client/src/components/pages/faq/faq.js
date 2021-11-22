import React from 'react';
import { Container, ListGroup, ListGroupItem } from 'react-bootstrap';
import { HashLink } from 'react-router-hash-link';

export default function Faq() {
  const sections = [
    {
      title: 'What public data is available in mSigPortal?',
      id: 'public-data',
      component: (
        <p>
          Collected public data on mSigPortal includes sample level mutational
          profiles, mutational signatures decomposition and associated variables
          collected from the scientific literature, such as{' '}
          <a href="https://dcc.icgc.org/pcawg" target="_blank">
            PCAWG
          </a>
          ,{' '}
          <a href="https://www.nature.com/articles/nature17676" target="_blank">
            Breast560
          </a>
          ,{' '}
          <a
            href="https://dceg.cancer.gov/research/cancer-types/lung/sherlock-lung-study"
            target="_blank"
          >
            Sherlock-Lung
          </a>
          ,{' '}
          <a
            href="https://www.cancer.gov/news-events/press-releases/2021/genetic-effects-chernobyl-radiation-exposure"
            target="_blank"
          >
            ChernobylThyroid
          </a>
          ,{' '}
          <a href="https://portal.gdc.cancer.gov/" target="_blank">
            TCGA
          </a>
          , and many others. Experimental strategies within these public data
          sources include mostly whole-genome sequencing (WGS) and whole-exome
          sequencing (WES) (currently only includes TCGA). We will continue to
          update the data as more become available.{' '}
        </p>
      ),
    },
    {
      title:
        'What reference signatures are included for use in analyses in mSigPortal?',
      id: 'reference-signatures',
      component: (
        <p>
          Several reference signature sets are included in mSigPortal. These
          include: Cancer Reference Signatures (SBS and RS),{' '}
          <a href="https://cancer.sanger.ac.uk/signatures/" target="_blank">
            COSMIC Signature
          </a>
          , (SBS, DBS, ID),{' '}
          <a
            href="https://signal.mutationalsignatures.com/explore/mutagens"
            target="_blank"
          >
            Environmental Mutagen Signatures
          </a>{' '}
          (SBS, ID),{' '}
          <a
            href="https://signal.mutationalsignatures.com/explore/study/1"
            target="_blank"
          >
            Organ-specific Cancer Signatures
          </a>{' '}
          (SBS and RS), SignatureAnalyzer PCAWG WGS Signatures (SBS, SBS, DBS,
          ID), SigProfiler PCAWG Strand Signatures (SBS), SigProfiler PCAWG WXS
          Signatures (SBS), and Other Published Signatures (SBS, ID). Click{' '}
          <a href="/#/catalog/signature">here</a> for additional information
          regarding current reference signatures in mSigPortal.
        </p>
      ),
    },
    {
      title: 'What are Single Base Substitutions (SBS)?',
      id: 'sbs',
      component: (
        <p>
          Single Base Substitutions (SBS) are when one base is substituted for
          another in a DNA sequence. In the context of mutational signatures,
          there are six classes of base substitution: C>A, C>G, C>T, T>A, T>C,
          and T>G. When also considering the bases to the immediate 5’ and 3’
          positions of the substitution, there are 96 possible mutations in this
          classification. When considering two bases 5’ and two bases 3’ of a
          mutation, the result is a 1536 channel matrix. For more information on
          SBS classification, click{' '}
          <a
            href="https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6041-2"
            target="_blank"
          >
            here
          </a>
          .
        </p>
      ),
    },
    {
      title: 'What are Doublet Base Substitutions (DBS)?',
      id: 'dbs',
      component: (
        <p>
          Doublet Base Substitutions (DBS) are when a set of two adjacent
          DNA-base pairs are substituted for another set of two adjacent DNA
          base pairs. There are 78 distinct DBS mutations categories, resulting
          in a 78 channel matrix. When considering the immediate 5’ base and 3’
          base of the DBS mutation, the matrix grows to 1248 channels. For more
          information on DBS classification, click{' '}
          <a
            href="https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6041-2"
            target="_blank"
          >
            here
          </a>
          .
        </p>
      ),
    },
    {
      title: 'What are Small Insertions and Deletions- Indels (ID)?',
      id: 'id',
      component: (
        <p>
          Small Insertions and Deletions (ID) are an incorporation of an
          additional set of base pairs, or a removal of a set of existing base
          pairs from a given location on a chromosome, respectively. They are
          typically {'<'}100 base pairs in length. In complicated events, the
          observed result is both a set of deleted base pairs and a set of
          inserted base pairs. Indels can be classified based on their length,
          nucleotides affected, and whether they occur at repetitive or
          microhomology regions. For more information on indels, click{' '}
          <a
            href="https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6041-2"
            target="_blank"
          >
            here
          </a>
          .
        </p>
      ),
    },
    {
      title:
        'What are homopolymers, repeat units, and microhomology for ID profile?',
      id: 'hrm',
      component: (
        <p>
          Homopolymers are repetitive regions that are stretches of the same
          base pair. Repeat units are repetitive regions of the same sequence of
          base pairs. Microhomology is when a long indel contains partially
          overlapping sequences. In other words, it is a degree of sequence
          similarity between the indel motif and the immediate junction
          sequence.
        </p>
      ),
    },
    {
      title: 'What is cosine similarity?',
      id: 'cosine-similarity',
      component: (
        <p>
          Cosine similarity is the similarity of signatures to one another. A
          cosine similarity value of 1 indicates that the signatures are exactly
          the same. A cosine similarity value of 0 indicates that the signatures
          are independent. Read more about the specifics of the formula for
          cosine similarity{' '}
          <a
            href="https://www.ncbi.nlm.nih.gov/research/mutagene/help#signatures"
            target="_blank"
          >
            here
          </a>
          .
        </p>
      ),
    },
    {
      title: 'What is Residual Sum of Squares RSS (RSS)?',
      id: 'rss',
      component: (
        <p>
          Residual Sum of Squares (RSS) is the measure of discrepancy between
          two profiles. Smaller RSS values indicate a greater similarity between
          signatures.
        </p>
      ),
    },
    {
      title: 'What is KL-Divergence?',
      id: 'kld',
      component: (
        <p>
          KL-Divergence is known as the Kullback-Leibler divergence, which
          measures the difference between two probability distributions for the
          same variable. The lower this value is, the better distributions
          match. For additional, detailed information, click{' '}
          <a
            href="https://www.ncbi.nlm.nih.gov/research/mutagene/help#signatures"
            target="_blank"
          >
            here
          </a>
          .
        </p>
      ),
    },
    {
      title: 'What are L1 norm and L2 norm?',
      id: 'l1-l2',
      component: (
        <p>
          L1 norm is the sum of the magnitudes of vectors in a space, and how
          one measures the distance between vectors. L2 norm is known as the
          Euclidean norm. It is the shortest distance to go from one point to
          another. For additional, detailed information, click{' '}
          <a
            href="https://www.kaggle.com/residentmario/l1-norms-versus-l2-norms"
            target="_blank"
          >
            here
          </a>
          .
        </p>
      ),
    },
    {
      title:
        'What is the meaning of IUPAC codes for nucleotide bases indicating specific mutational types?',
      id: 'iupac',
      component: (
        <div>
          <p>
            Traditionally, we see nucleotide bases in the notation of A, T, C,
            and G. However, there are often multiple possibilities as to the
            base to the immediate 5’ and 3’ position of the mutated base. The
            International Union of Pure and Applied Chemistry (IUPAC) has a
            system for giving codes to identify the nucleotide bases. Use the
            following table to investigate the different possibilities and the
            letters they are represented by:
          </p>
          <p>
            Example: Using the table below, a mutational pattern represented by
            HCN>HTN is a C>T mutation. The H indicates that the base in this
            position could be A, C, or T. The N indicates that the base in this
            position could be A, T, C, or G.
          </p>
          <table className="w-auto table table-bordered table-striped text-center">
            <thead>
              <tr>
                <th>
                  Letter to Represent
                  <br />
                  the Base Possibilities
                </th>
                <th>
                  Base
                  <br />
                  Possibilities
                </th>
              </tr>
            </thead>
            <tbody>
              <tr>
                <td>R</td>
                <td>A or G</td>
              </tr>
              <tr>
                <td>Y</td>
                <td>C or T</td>
              </tr>
              <tr>
                <td>S</td>
                <td>G or C</td>
              </tr>
              <tr>
                <td>W</td>
                <td>A or T</td>
              </tr>
              <tr>
                <td>K</td>
                <td>G or T</td>
              </tr>
              <tr>
                <td>M</td>
                <td>A or C</td>
              </tr>
              <tr>
                <td>B</td>
                <td>C, G, or T</td>
              </tr>
              <tr>
                <td>D</td>
                <td>A, G, or T</td>
              </tr>
              <tr>
                <td>H</td>
                <td>A, C, or T</td>
              </tr>
              <tr>
                <td>V</td>
                <td>A, C, or G</td>
              </tr>
              <tr>
                <td>N</td>
                <td>A, T, C, or G</td>
              </tr>
            </tbody>
          </table>
        </div>
      ),
    },
    {
      title: 'What is mutational pattern enrichment analysis?',
      id: 'mpea',
      component: (
        <p>
          This type of analysis aims to determine frequency and enrichment of
          different mutational patterns considering the IUPAC codes of the bases
          to the immediate 5’ and 3’ positions of the substitution. For
          enrichment analysis, mSigPortal calculates and visualizes the
          proportion of specific mutation pattern context (such as HCG>HTG) as
          compared to other contexts with the same mutation type (such as C>G).
          This analysis will help to identify the samples with specific
          mutational context due to the same mutational process.
        </p>
      ),
    },
    {
      title: 'What is Principal Component Analysis (PCA)?',
      id: 'pca',
      component: (
        <p>
          Principal Component Analysis (PCA) is a dimensionality-reduction
          method that helps to explain the variation found in the data through
          the establishment of different principal components. Principal
          components are new variables that are a linear combination of the
          initial variables. These are uncorrelated to one another, and the
          attempt is to obtain the maximum amount of variance within the first
          components. The number of dimensions in the data equals the number of
          principal components generated. In the case of mutational signatures,
          each principal component can also be used to compare known mutational
          signatures{' '}
          <a href="https://www.nature.com/articles/nbt0308-303" target="_blank">
            (Ringnér, 2008)
          </a>
          .
        </p>
      ),
    },
    {
      title: 'What is kataegis?',
      id: 'kataegis',
      component: (
        <p>
          Kataegis is localized substitution hypermutation, often characterized
          by clusters of C>T and/or C>G mutations, commonly at TpCpN
          trinucleotides. While it is common to think that mutations would occur
          every million base pairs or so, these localized instances of
          substitution hypermutations occur at much shorter base pair distances.
          This occurrence can be visualized on a rainfall plot, leading to
          kataegis getting its name from the greek word for thunderstorm.
        </p>
      ),
    },
    {
      title: 'What is tumor mutational burden (TMB)?',
      id: 'tmb',
      component: (
        <p>
          Tumor mutational burden (TMB) is the number of mutations in the DNA of
          cancer cells. TMB is calculated by the number of non-synonymous
          somatic mutations (single nucleotide variants and small
          insertions/deletions) per mega-base{' '}
          <a
            href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5461196/"
            target="_blank"
          >
            (Zehir et al., 2017)
          </a>
          . Tumors with high mutational burden may respond to treatments that
          fall into the category of immunotherapy{' '}
          <a
            href="https://www.cancer.gov/publications/dictionaries/cancer-terms/def/tumor-mutational-burden"
            target="_blank"
          >
            (NCI, 2021)
          </a>
          .
        </p>
      ),
    },
    {
      title: 'What is mutational signature burden (MSB)?',
      id: 'msb',
      component: (
        <p>
          Mutational signature burden (MSB) is the number of mutations
          associated with a given mutational signature in the DNA of cancer
          cells.
        </p>
      ),
    },
    {
      title:
        'What types of association tests can be performed using mSigPortal?',
      id: 'association',
      component: (
        <p>
          Currently, only univariable and multivariable association tests can be
          performed using mSigPortal. Univariable association tests investigate
          relationships between a variable from provided data and a signature
          exposure variable using different statistical tests depending on the
          variables type. Multivariable association tests perform linear or
          logistic regression based on data from the selected variables and a
          signature exposure variable.
        </p>
      ),
    },
  ];

  return (
    <div className="mx-3">
      <div className="bg-white border p-3 mx-3">
        <Container>
          <div id="help-logo" className="text-center">
            <img
              height="150"
              src="assets/images/msigportal-logo.png"
              alt="mSigPortal Logo"
            />
          </div>
          <div className="mt-3">
            <p>
              mSigPortal is designed to be an easy and intuitive web portal to
              explore, visualize and analyze mutational signature related data
              for genomic studies. This documentation page collects the most
              frequently asked questions about mSigPortal.
            </p>
            <ListGroup>
              {sections.map(({ title, id }) => (
                <ListGroupItem>
                  <HashLink
                    smooth
                    to={{
                      pathname: '/faq',
                      hash: `#${id}`,
                      state: { fromDashboard: true },
                    }}
                  >
                    {title}
                  </HashLink>
                </ListGroupItem>
              ))}
            </ListGroup>
          </div>
          {sections.map(({ title, id, component }) => (
            <div id={id} className="my-5">
              <p>
                <b>{title}</b>
              </p>
              {component}
            </div>
          ))}
        </Container>
      </div>
    </div>
  );
}
