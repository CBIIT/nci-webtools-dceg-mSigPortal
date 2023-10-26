export const schema = [
  {
    name: 'association',
    schema: (table) => {
      table.increments('id');
      table.string('study');
      table.string('strategy');
      table.string('cancer');
      table.string('sample');
      table.string('icgcSpecimenId');
      table.string('icgcDonorId');
      table.string('dataSource');
      table.string('dataType');
      table.string('variableName');
      table.string('variableValue');
      table.string('variableValueType');
    },
    index: (table) => {
      table.index(['study', 'strategy', 'cancer', 'sample']);
    },
  },

  {
    name: 'associationOption',
    type: 'materializedView',
    dependsOn: ['association'],
    schema: (view, connection) => {
      const columns = ['study', 'strategy', 'cancer', 'sample'];
      const query = connection('association').distinct(columns);
      view.columns(columns);
      return view.as(query);
    },
    index: (table) => {
      table.index(['study', 'strategy', 'cancer', 'sample']);
    },
  },

  {
    name: 'exposure',
    schema: (table) => {
      table.increments('id');
      table.string('study');
      table.string('strategy');
      table.string('cancer');
      table.string('organ');
      table.string('sample');
      table.string('signatureSetName');
      table.string('signatureName');
      table.double('exposure');
    },
    index: (table) => {
      table.index(['study', 'strategy', 'cancer', 'signatureSetName']);
    },
  },

  {
    name: 'exposureOption',
    type: 'materializedView',
    dependsOn: ['exposure'],
    schema: (view, connection) => {
      const columns = [
        'study',
        'strategy',
        'cancer',
        'signatureSetName',
        'signatureName',
      ];
      const query = connection('exposure').distinct(columns);
      view.columns(columns);
      view.as(query);
    },
    index: (table) => {
      table.index([
        'study',
        'strategy',
        'cancer',
        'signatureSetName',
        'signatureName',
      ]);
    },
  },

  {
    name: 'seqmatrix',
    schema: (table) => {
      table.increments('id');
      table.string('study');
      table.string('strategy');
      table.string('cancer');
      table.string('sample');
      table.string('profile');
      table.integer('matrix');
      table.string('mutationType');
      table.integer('mutations');
    },
    index: (table) => {
      table.index([
        'study',
        'strategy',
        'cancer',
        'sample',
        'profile',
        'matrix',
      ]);
    },
  },

  {
    name: 'seqmatrixOption',
    type: 'materializedView',
    dependsOn: ['seqmatrix'],
    schema: (view, connection) => {
      const columns = [
        'study',
        'strategy',
        'cancer',
        'sample',
        'profile',
        'matrix',
      ];
      const query = connection('seqmatrix').distinct(columns);
      view.columns(columns);
      return view.as(query);
    },
    index: (table) => {
      table.index([
        'study',
        'strategy',
        'cancer',
        'sample',
        'profile',
        'matrix',
      ]);
    },
  },

  {
    name: 'seqmatrixSummary',
    type: 'materializedView',
    dependsOn: ['seqmatrix'],
    schema: (view, connection) => {
      const columns = [
        'study',
        'strategy',
        'cancer',
        'sample',
        'profile',
        'matrix',
        'logTotalMutations',
        'meanTotalMutations',
      ];
      const totalCountColumns = [
        ...columns.slice(0, -2),
        connection.raw('sum(mutations) as "totalMutations"'),
      ];
      const summaryColumns = [
        ...columns.slice(0, -3),
        connection.raw(
          `string_agg(cast(matrix as text), '/' order by matrix asc) as matrix`
        ),
        connection.raw(`log10("totalMutations") as "logTotalMutations"`),
        connection.raw(`avg("totalMutations") as "meanTotalMutations"`),
      ];
      const totalCountQuery = connection
        .select(totalCountColumns)
        .from('seqmatrix')
        .groupBy(totalCountColumns.slice(0, -1))
        .orderBy(totalCountColumns.slice(0, -1).map((column) => ({ column })));
      const summaryQuery = connection
        .with('totalCounts', totalCountQuery)
        .select(summaryColumns)
        .from('totalCounts')
        .where('totalMutations', '>', 0)
        .groupBy([...summaryColumns.slice(0, -3), 'totalMutations']);

      view.columns(columns);
      return view.as(summaryQuery);
    },
    index: (table) => {
      table.index([
        'study',
        'strategy',
        'cancer',
        'sample',
        'profile',
        'matrix',
      ]);
    },
  },

  {
    name: 'signature',
    schema: (table) => {
      table.increments('id');
      table.string('source');
      table.string('strategy');
      table.string('profile');
      table.integer('matrix');
      table.string('signatureSetName');
      table.string('strandInfo');
      table.string('strand');
      table.string('signatureName');
      table.string('mutationType');
      table.double('contribution');
    },
    index: (table) => {
      table.index([
        'source',
        'strategy',
        'profile',
        'matrix',
        'signatureSetName',
        'signatureName',
      ]);
    },
  },

  {
    name: 'cache',
    schema: (table) => {
      table.text('key');
      table.json('value');
    },
    index: (table) => {
      table.unique([
        'key',
      ]);
    },
  },

  {
    name: 'signatureOption',
    type: 'materializedView',
    dependsOn: ['signature'],
    schema: (view, connection) => {
      const columns = [
        'source',
        'strategy',
        'profile',
        'matrix',
        'signatureSetName',
        'signatureName',
      ];
      const query = connection('signature').distinct(columns);
      view.columns(columns);
      return view.as(query);
    },
    index: (table) => {
      table.index([
        'source',
        'strategy',
        'profile',
        'matrix',
        'signatureSetName',
        'signatureName',
      ]);

      table.index(['signatureSetName', 'signatureName']);
    },
  },

  {
    name: 'signatureSummary',
    schema: (table) => {
      table.increments('id');
      table.string('species');
      table.string('profile');
      table.integer('matrix');
      table.string('signatureSetName');
      table.string('count');
    },
    index: (table) => {
      table.index(['profile', 'matrix', 'signatureSetName']);
    },
  },

  {
    name: 'etiologyOptions',
    schema: (table) => {
      table.increments('id');
      table.string('category');
      table.string('etiology');
      table.string('signature');
      table.json('json');
    },
    index: (table) => {
      table.index(['category']);
    },
  },

  {
    name: 'etiology',
    schema: (table) => {
      table.increments('id');
      table.string('study');
      table.string('strategy');
      table.string('cancer');
      table.string('organ');
      table.string('sample');
      table.string('signatureSetName');
      table.integer('mutations');
      table.double('cosineSimilarity');
      table.integer('sampleSize');
      table.string('signatureName');
      table.double('exposure');
      table.double('burden');
      table.integer('signatureSize');
    },
    index: (table) => {
      table.index(['study', 'strategy', 'signatureName', 'signatureSetName']);
    },
  },

  {
    name: 'etiologyOrgan',
    schema: (table) => {
      table.increments('id');
      table.string('signature');
      table.string('cohort');
      table.string('organ');
      table.string('prevalence');
      table.string('organSpecificSignature');
      table.string('contribution');
    },
    index: (table) => {
      table.index(['signature', 'cohort', 'organ']);
    },
  },

  {
    name: 'publication',
    schema: (table) => {
      table.increments('id');
      table.string('category');
      table.string('firstAuthor');
      table.string('lastAuthor');
      table.integer('year');
      table.string('journal');
      table.string('bioRxivOrPubmedId');
      table.text('title');
      table.string('doi');
      table.string('note');
      table.string('diseaseOrPhenotypeOrExposure');
      table.string('cancerType');
      table.string('experimentalStrategy');
      table.string('softwareName');
      table.string('computationalMethod');
      table.string('programmingLanguage');
      table.string('sourceUrl');
    },
    index: (table) => {
      table.index(['category']);
    },
  },

  {
    name: 'pattern',
    schema: (table) => {
      table.increments('id');
      table.string('study');
      table.string('cancer');
      table.string('sample');
      table.integer('total');
      table.string('pattern');
      table.integer('n0');
      table.double('n1');
      table.double('n2');
    },
    index: (table) => {
      table.index(['study', 'cancer', 'n1']);
    },
  },

  {
    name: 'refgenome',
    schema: (table) => {
      table.increments('id');
      table.string('genome');
      table.string('chr');
      table.double('len');
      table.double('start');
      table.double('end');
    },
    index: (table) => {
      table.index(['genome']);
    },
  },

  {
    name: 'importLog',
    schema: (table, connection) => {
      table.increments('id');
      table.string('status');
      table.text('log');
      table.timestamp('createdAt').defaultTo(connection.fn.now());
      table.timestamp('updatedAt').defaultTo(connection.fn.now());
    },
  },
];
