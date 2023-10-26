const schema = [
  {
    name: 'seqmatrix',
    schema: (table) => {
      table.increments('id');
      table.string('sample');
      table.string('mutationType');
      table.integer('mutations');
      table.string('cancer');
    },
    index: (table) => {
      table.index(['sample']);
    },
  },
  {
    name: 'signature',
    schema: (table) => {
      table.increments('id');
      table.string('signatureName');
      table.string('mutationType');
      table.double('contribution');
    },
    index: (table) => {
      table.index(['signatureName']);
    },
  },
  {
    name: 'exposure',
    schema: (table) => {
      table.increments('id');
      table.string('sample');
      table.string('signatureName');
      table.double('exposure');
      table.string('cancer');
    },
    index: (table) => {
      table.index(['sample']);
    },
  },
  {
    name: 'seqmatrixOption',
    type: 'materializedView',
    dependsOn: ['seqmatrix'],
    create: (connection) => {
      const columns = ['sample'];
      return connection.raw(
        [
          `CREATE TABLE seqmatrixOption AS`,
          `SELECT DISTINCT ${columns.join(',')}`,
          `FROM seqmatrix`,
        ].join(' ')
      );
    },
    index: (table) => {
      table.index(['sample']);
    },
  },
  {
    name: 'signatureOption',
    type: 'materializedView',
    dependsOn: ['signature'],
    create: (connection) => {
      const columns = ['signatureName'];
      return connection.raw(
        [
          `CREATE TABLE signatureOption AS`,
          `SELECT DISTINCT ${columns.join(',')}`,
          `FROM signature`,
        ].join(' ')
      );
    },
    index: (table) => {
      table.index(['signatureName']);
    },
  },
  {
    name: 'exposureOption',
    type: 'materializedView',
    dependsOn: ['exposure'],
    create: (connection) => {
      const columns = ['sample'];
      return connection.raw(
        [
          `CREATE TABLE exposureOption AS`,
          `SELECT DISTINCT ${columns.join(',')}`,
          `FROM exposure`,
        ].join(' ')
      );
    },
    index: (table) => {
      table.index(['sample']);
    },
  },
];
export { schema };
