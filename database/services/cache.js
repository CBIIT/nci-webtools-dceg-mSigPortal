// The server caches /treeLeaf responses (built from signature, exposure and seqmatrix) with no expiry,
// and imports never recreate the cache table, so stale results would outlive the data they came from.
export const clearResponseCache = {
  description: "Clear cached API responses",
  type: "postImport",
  callback: async (connection, logger) => {
    const { rows } = await connection.query(
      `select to_regclass('"cache"') is not null as "exists"`,
    );
    if (!rows[0]?.exists) {
      logger.info(`No cache table found, skipping cache clear`);
      return;
    }
    await connection.query(`truncate "cache"`);
    logger.info(`Cleared cached API responses`);
  },
};
