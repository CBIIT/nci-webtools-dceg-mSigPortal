// Reference signatures have no study in their CSV; client queries filter on study = 'Reference'.
// Must run before the signatureOption view is refreshed.
export const backfillReferenceSignatureStudy = {
  description: "Backfill study for reference signatures",
  type: "postImport",
  callback: async (connection, logger) => {
    const start = Date.now();
    try {
      // de novo rows already carry their own study and must not be relabelled
      const { rowCount } = await connection.query(
        `update "signature" set "study" = 'Reference' where "study" is null and "source" <> 'Study_signatures'`,
      );
      const duration = ((Date.now() - start) / 1000).toFixed(2);
      logger.info(
        `Backfilled study for ${rowCount} reference signature rows in ${duration}s`,
      );
    } catch (error) {
      const duration = ((Date.now() - start) / 1000).toFixed(2);
      logger.error(
        `Failed to backfill study for reference signatures after ${duration}s: ${error.message}`,
      );
      throw error;
    }
  },
};
