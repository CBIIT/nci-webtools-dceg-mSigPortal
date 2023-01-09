export const sources = [
  {
    sourcePath: 'Association/data.csv',
    table: 'association',
    description: 'association data',
    columns: [
      'study',
      'strategy',
      'cancer',
      'sample',
      'icgcSpecimenId',
      'icgcDonorId',
      'dataSource',
      'dataType',
      'variableName',
      'variableValue',
      'variableValueType',
    ],
  },

  {
    sourcePath: 'Exposure/data.csv',
    table: 'exposure',
    description: 'exposure data',
    columns: [
      'study',
      'strategy',
      'cancer',
      'organ',
      'sample',
      'signatureSetName',
      'signatureName',
      'exposure',
    ],
  },

  {
    sourcePath: 'Seqmatrix/data.csv',
    table: 'seqmatrix',
    description: 'seqmatrix data',
    columns: [
      'study',
      'cancer',
      'sample',
      'strategy',
      'profile',
      'matrix',
      'mutationType',
      'mutations',
    ],
  },

  {
    sourcePath: 'Signature/data.csv',
    table: 'signature',
    description: 'signature data',
    columns: [
      'source',
      'profile',
      'matrix',
      'signatureSetName',
      'strategy',
      'strandInfo',
      'strand',
      'signatureName',
      'mutationType',
      'contribution',
    ],
  },

  {
    sourcePath: 'Signature/summary.csv',
    table: 'signatureSummary',
    description: 'signature summary',
    columns: ['species', 'profile', 'matrix', 'signatureSetName', 'count'],
  },

  {
    sourcePath: 'Etiology/Etiology_cancer_specific_signatures_2022.csv',
    table: 'etiologyOptions',
    description: 'etiology - cancer_specific_signatures_2022',
    columns: ['category', 'etiology', 'signature', 'json'],
  },

  {
    sourcePath: 'Etiology/Etiology_cancer_specific_signatures.csv',
    table: 'etiologyOptions',
    description: 'etiology - cancer_specific_signatures',
    columns: ['category', 'etiology', 'signature', 'json'],
  },

  {
    sourcePath: 'Etiology/Etiology_cancer_therapies.csv',
    table: 'etiologyOptions',
    description: 'etiology - cancer_therapies',
    columns: ['category', 'etiology', 'signature', 'json'],
  },
  {
    sourcePath: 'Etiology/Etiology_cosmic.csv',
    table: 'etiologyOptions',
    description: 'etiology - cosmic',
    columns: ['category', 'etiology', 'signature', 'json'],
  },
  {
    sourcePath: 'Etiology/Etiology_enviromental_mutagenesis.csv',
    table: 'etiologyOptions',
    description: 'etiology - enviromental_mutagenesis',
    columns: ['category', 'etiology', 'signature', 'json'],
  },
  {
    sourcePath: 'Etiology/Etiology_gene_edits.csv',
    table: 'etiologyOptions',
    description: 'etiology - gene_edits',
    columns: ['category', 'etiology', 'signature', 'json'],
  },

  {
    sourcePath: 'Etiology/Etiology_others.csv',
    table: 'etiologyOptions',
    description: 'etiology - others',
    columns: ['category', 'etiology', 'signature', 'json'],
  },

  {
    sourcePath: 'Etiology/etiology.csv',
    table: 'etiology',
    description: 'etiology - aetiology_exposure',
    columns: [
      'study',
      'strategy',
      'cancer',
      'organ',
      'sample',
      'signatureSetName',
      'mutations',
      'cosineSimilarity',
      'sampleSize',
      'signatureName',
      'exposure',
      'burden',
      'signatureSize',
    ],
  },

  {
    sourcePath: 'Etiology/etiologyOrgan.csv',
    table: 'etiologyOrgan',
    description: 'etiologyOrgan - aetiology_organ_specific_signature',
    columns: [
      'signature',
      'cohort',
      'organ',
      'prevalence',
      'organSpecificSignature',
      'contribution',
    ],
  },

  {
    sourcePath: 'Others/Publications.csv',
    table: 'publication',
    description: 'publications',
    columns: [
      'category',
      'diseaseOrPhenotypeOrExposure',
      'cancerType',
      'experimentalStrategy',
      'firstAuthor',
      'lastAuthor',
      'year',
      'journal',
      'bioRxivOrPubmedId',
      'title',
      'doi',
      'note',
      'softwareName',
      'computationalMethod',
      'programmingLanguage',
      'sourceUrl',
    ],
  },

  {
    sourcePath: 'Others/pattern.csv',
    table: 'pattern',
    description: 'pattern',
    columns: [
      'study',
      'cancer',
      'sample',
      'total',
      'pattern',
      'n0',
      'n1',
      'n2',
    ],
  },

  {
    sourcePath: 'Others/refgenome.csv',
    table: 'refgenome',
    description: 'genome - chromosome sizes',
    columns: ['genome', 'chr', 'len', 'start', 'end'],
  },

  {
    description: 'Refresh materialized views',
    type: 'postImport',
    callback: async (connection) => {
      await connection.query('refresh materialized view "exposureOption"');
      await connection.query('refresh materialized view "seqmatrixOption"');
      await connection.query('refresh materialized view "seqmatrixSummary"');
      await connection.query('refresh materialized view "signatureOption"');
      await connection.query('refresh materialized view "associationOption"');
    },
  },
];
