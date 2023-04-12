// imports json data into sqlite file
export function sqliteImport(connection, data, userSchema) {
  return new Promise(async (resolve, reject) => {
    const tables = userSchema.filter((e) => !e.type || e.type === 'table');
    const materializedViews = userSchema.filter(
      (e) => e.type === 'materializedView'
    );
    const indexedTables = userSchema.filter(
      (s) => typeof s.index === 'function'
    );
    try {
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
      resolve(true);
    } catch (error) {
      reject(error);
    }
  });
}
