

function createDatabaseCache(connection, tableName = 'cache')  {
  return {
    get: async (key) => {
      const result = await connection(tableName).where({ key }).first();
      return result ? result.value : null;
    },
    set: async (key, value) => {
      await connection(tableName).insert({ key, value }).onConflict('key').merge();
    },
    clear: async (key) => {
      await connection(tableName).where({ key }).del();
    },
    reset: async () => {
      await connection(tableName).truncate();
    }
  }
}
