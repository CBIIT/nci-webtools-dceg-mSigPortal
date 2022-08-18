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
    schema: (view, connection) => {
      const columns = ['profile', 'matrix', 'sample'];
      const query = connection('seqmatrix').distinct(columns);
      return view.as(query);
    },
    index: (table) => {
      table.index(['profile', 'matrix', 'sample']);
    },
  },

  //   {
  //     name: 'seqmatrixSummary',
  //     type: 'materializedView',
  //     dependsOn: ['seqmatrix'],
  //     create: (connection) => {
  //       const columns = [
  //         'study',
  //         'strategy',
  //         'cancer',
  //         'sample',
  //         'profile',
  //         'matrix',
  //         'logTotalMutations',
  //         'meanTotalMutations',
  //       ];

  //       const totalCountColumns = [
  //         ...columns.slice(0, -2),
  //         connection.raw('sum(mutations) as "totalMutations"'),
  //       ];
  //       const summaryColumns = [
  //         ...columns.slice(0, -3),
  //         connection.raw(
  //           `string_agg(cast(matrix as text), '/' order by matrix asc) as matrix`
  //         ),
  //         connection.raw(`log10(sum("totalMutations")) as "logTotalMutations"`),
  //         connection.raw(`avg("totalMutations") as "meanTotalMutations"`),
  //       ];
  //       const totalCountQuery = connection
  //         .select(totalCountColumns)
  //         .from('seqmatrix')
  //         .groupBy(totalCountColumns.slice(0, -1))
  //         .orderBy(totalCountColumns.slice(0, -1).map((column) => ({ column })));
  //       const summaryQuery = connection
  //         .with('totalCounts', totalCountQuery)
  //         .select(summaryColumns)
  //         .from('totalCounts')
  //         .where('totalMutations', '>', 0)
  //         .groupBy([...summaryColumns.slice(0, -3), 'totalMutations']);

  //       view.columns(columns);
  //       return view.as(summaryQuery);
  //     },
  //     schema: (view, connection) => {
  //       const columns = [
  //         'study',
  //         'strategy',
  //         'cancer',
  //         'sample',
  //         'profile',
  //         'matrix',
  //         'logTotalMutations',
  //         'meanTotalMutations',
  //       ];
  //       const totalCountColumns = [
  //         ...columns.slice(0, -2),
  //         connection.raw('sum(mutations) as "totalMutations"'),
  //       ];
  //       const summaryColumns = [
  //         ...columns.slice(0, -3),
  //         connection.raw(
  //           `string_agg(cast(matrix as text), '/' order by matrix asc) as matrix`
  //         ),
  //         connection.raw(`log10(sum("totalMutations")) as "logTotalMutations"`),
  //         connection.raw(`avg("totalMutations") as "meanTotalMutations"`),
  //       ];
  //       const totalCountQuery = connection
  //         .select(totalCountColumns)
  //         .from('seqmatrix')
  //         .groupBy(totalCountColumns.slice(0, -1))
  //         .orderBy(totalCountColumns.slice(0, -1).map((column) => ({ column })));
  //       const summaryQuery = connection
  //         .with('totalCounts', totalCountQuery)
  //         .select(summaryColumns)
  //         .from('totalCounts')
  //         .where('totalMutations', '>', 0)
  //         .groupBy([...summaryColumns.slice(0, -3), 'totalMutations']);

  //       view.columns(columns);
  //       return view.as(summaryQuery);
  //     },
  //     index: (table) => {
  //       table.index([
  //         'study',
  //         'strategy',
  //         'cancer',
  //         'sample',
  //         'profile',
  //         'matrix',
  //       ]);
  //     },
  //   },
];

module.exports = { schema };
