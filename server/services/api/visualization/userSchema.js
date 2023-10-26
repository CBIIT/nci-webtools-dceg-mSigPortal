const schema = [
  {
    name: 'seqmatrix',
    schema: (table) => {
      table.increments('id');
      table.string('sample');
      table.string('profile');
      table.integer('matrix');
      table.string('mutationType');
      table.integer('mutations');
      table.string('filter');
    },
    index: (table) => {
      table.index(['profile', 'matrix', 'sample', 'filter']);
    },
  },
  {
    name: 'cluster',
    schema: (table) => {
      table.increments('id');
      table.string('sample');
      table.string('geneId');
      table.string('genome');
      table.string('mutType');
      table.string('chr');
      table.integer('start');
      table.integer('end');
      table.string('ref');
      table.string('alt');
      table.string('mutClass');
      table.integer('IMDplot');
      table.integer('group');
      table.integer('IMD');
      table.double('VAF/CCF');
      table.string('subclass');
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
      const columns = ['profile', 'matrix', 'sample', 'filter'];
      return connection.raw(
        [
          `CREATE TABLE seqmatrixOption AS`,
          `SELECT DISTINCT ${columns.join(',')}`,
          `FROM seqmatrix`,
        ].join(' ')
      );
    },
    index: (table) => {
      table.index(['profile', 'matrix', 'sample', 'filter']);
    },
  },
  {
    name: 'seqmatrixSummary',
    type: 'materializedView',
    dependsOn: ['seqmatrix'],
    create: (connection) => {
      return connection.raw(`
          CREATE TABLE seqmatrixSummary AS
          WITH records AS (
            SELECT        
              sample || IIF(filter != '', '@', '') || COALESCE(filter, '') AS sample,
              profile,
              matrix,
              filter,
              sum(mutations) as totalMutations
            FROM seqmatrix
            GROUP BY sample, filter, profile, matrix
          )
          SELECT
            sample,
            profile,
            filter,
            GROUP_CONCAT(cast(matrix as text), '/') AS matrix,
            log10(totalMutations) as logTotalMutations,
            avg(totalMutations) as meanTotalMutations
          FROM records r
          WHERE totalMutations > 0
          GROUP BY sample, profile, totalMutations
          ORDER BY profile
          `);
    },
    index: (table) => {
      table.index(['sample', 'profile', 'matrix', 'filter']);
    },
  },
];

export { schema };
