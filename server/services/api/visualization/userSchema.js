const schema = [
  {
    name: 'seqmatrix',
    schema: (table) => {
      table.increments('id');
      table.string('profile');
      table.integer('matrix');
      table.string('sample');
      table.string('mutationType');
      table.integer('mutations');
    },
    index: (table) => {
      table.index(['profile', 'matrix', 'sample']);
    },
  },

  {
    name: 'seqmatrixOption',
    type: 'materializedView',
    dependsOn: ['seqmatrix'],
    create: (connection) => {
      const columns = ['profile', 'matrix', 'sample'];

      return connection.raw(
        [
          `CREATE TABLE seqmatrixOption AS`,
          `SELECT DISTINCT ${columns.join(',')}`,
          `FROM seqmatrix`,
        ].join(' ')
      );
    },
    index: (table) => {
      table.index(['profile', 'matrix', 'sample']);
    },
  },

  {
    name: 'seqmatrixSummary',
    type: 'materializedView',
    dependsOn: ['seqmatrix'],
    create: (connection) => {
      return connection.raw(
        `
          CREATE TABLE seqmatrixSummary AS
          WITH records AS (
            SELECT        
                sample,
                profile,
                matrix,
                sum(mutations) as totalMutations
            FROM seqmatrix
            GROUP BY sample, profile, matrix
          )
          SELECT
            sample,
            profile,
            group_concat(cast(matrix as text), '/') as matrix,
            log10(totalMutations) as logTotalMutations,
            avg(totalMutations) as meanTotalMutations
          FROM records r
          WHERE totalMutations > 0
          GROUP BY sample, profile, totalMutations
          ORDER BY profile
          `
      );
    },
    index: (table) => {
      table.index(['sample', 'profile', 'matrix']);
    },
  },
];

module.exports = { schema };
