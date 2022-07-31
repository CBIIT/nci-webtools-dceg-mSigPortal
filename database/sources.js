export const sources = [
  {
    sourcePath: "Association/data.csv",
    table: "association",
    description: "association data",
    columns: [
      "study",
      "strategy",
      "cancer",
      "sample",
      "icgcSpecimenId",
      "icgcDonorId",
      "dataSource",
      "dataType",
      "variableName",
      "variableValue",
      "variableValueType",
    ],
  },

  {
    sourcePath: "Exposure/data.csv",
    table: "exposure",
    description: "exposure data",
    columns: [
      "study",
      "strategy",
      "cancer",
      "organ",
      "sample",
      "signatureSetName",
      "signatureName",
      "exposure",
    ],
  },

  {
    sourcePath: "Seqmatrix/data.csv",
    table: "seqmatrix",
    description: "seqmatrix data",
    columns: [
      "study",
      "cancer",
      "sample",
      "strategy",
      "profile",
      "matrix",
      "mutationType",
      "mutations",
    ],
  },

  {
    sourcePath: "Signature/data.csv",
    table: "signature",
    description: "signature data",
    columns: [
      "source",
      "profile",
      "matrix",
      "signatureSetName",
      "strategy",
      "strandInfo",
      "strand",
      "signatureName",
      "mutationType",
      "contribution",
    ],
  },

  {
    sourcePath: "Etiology/Etiology_cancer_specific_signatures.csv",
    table: "etiology",
    description: "etiology - cancer_specific_signatures",
    columns: [
      "etiology",
      "tissueSpecificSignature",
      "refSignature",
      "refSignatureProportion",
      "study",
      "studyUrl",
      "source",
      "sourceUrl",
      "description",
      "category",
    ]
  },
  {
    sourcePath: "Etiology/Etiology_cancer_therapies.csv",
    table: "etiology",
    description: "etiology - cancer_therapies",
    columns: [
      "treatment",
      "signature",
      "signatureExtractionMethod",
      "tumorType",
      "study",
      "studyUrl",
      "category",
    ]
  },
  {
    sourcePath: "Etiology/Etiology_cosmic.csv",
    table: "etiology",
    description: "etiology - cosmic",
    columns: [
      "etiology",
      "signatureName",
      "signatureSource",
      "url",
      "description",
      "descriptionStrandBias",
      "study",
      "category",
    ]
  },
  {
    sourcePath: "Etiology/Etiology_enviromental_mutagenesis.csv",
    table: "etiology",
    description: "etiology - enviromental_mutagenesis",
    columns: [
      "etiology",
      "signature",
      "mutagen",
      "treatment",
      "study",
      "studyUrl",
      "source",
      "sourceUrl",
      "category",
    ]
  },
  {
    sourcePath: "Etiology/Etiology_gene_edits.csv",
    table: "etiology",
    description: "etiology - gene_edits",
    columns: [
      "etiology",
      "signature",
      "cellLine",
      "study",
      "studyUrl",
      "source",
      "sourceUrl",
      "category",
    ]
  },
 
  {
    sourcePath: "Etiology/Etiology_others.csv",
    table: "etiology",
    description: "etiology - others",
    columns: [
      "etiology",
      "signatureName",
      "signatureSource",
      "url",
      "description", 
      "category",
    ]
  },

  {
    sourcePath: "Others/Publications.csv",
    table: "publication",
    description: "publications",
    columns: [
      "category",
      "diseaseOrPhenotypeOrExposure",
      "cancerType",
      "experimentalStrategy",
      "firstAuthor",
      "lastAuthor",
      "year",
      "journal",
      "bioRxivOrPubmedId",
      "title",
      "doi",
      "note",
      "softwareName",
      "computationalMethod",
      "programmingLanguage",
      "sourceUrl",
    ]
  },

  {
    description: "Refresh materialized views",
    type: "postImport",
    callback: async (connection) => {
      await connection.query('REFRESH MATERIALIZED VIEW "seqmatrixOption"');
      await connection.query('REFRESH MATERIALIZED VIEW "seqmatrixSummary"');
    }
  }
];
