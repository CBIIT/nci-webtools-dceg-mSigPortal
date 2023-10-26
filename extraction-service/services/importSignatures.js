export const signatureSchema = [
  {
    name: 'signature',
    schema: (table) => {
      table.increments('id');
      table.string('signatureName');
      table.string('mutationType');
      table.integer('mutations');
    },
    index: (table) => {
      table.index(['signatureName']);
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
];

export async function importUserSession(
  connection,
  data,
  schema = signatureSchema
) {
  const tables = signatureSchema.filter((e) => !e.type || e.type === 'table');
  const materializedViews = schema.filter((e) => e.type === 'materializedView');
  const indexedTables = schema.filter((s) => typeof s.index === 'function');

  // create tables
  for (const { name, schema } of tables) {
    await connection.schema.createTable(name, (table) =>
      schema(table, connection)
    );
  }

  // import data
  for (const [tableName, tableData] of Object.entries(data))
    await connection.batchInsert(tableName, tableData, 100);
  // create "materialized" style tables

  for (const { create } of materializedViews) {
    await create(connection);
  }

  // index tables
  for (const { name, index } of indexedTables) {
    await connection.schema.table(name, index);
  }
  return true;
}
