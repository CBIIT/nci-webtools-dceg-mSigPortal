export const sources = [
  {
    sourcePath: "Signature/data.csv",
    table: "signature",
    description: "signature data",
    columns: [
      "source",
      "profile",
      "matrix",
      "signatureSetName",
      "strategy",
      "strandInfo",
      "strand",
      "signatureName",
      "mutationType",
      "contribution",
    ],
  },

  {
    sourcePath: "Exposure/Study_Signatures/data.csv",
    table: "signature",
    description: "study de novo signature data",
    columns: [
      "source",
      "profile",
      "matrix",
      "signatureSetName",
      "strategy",
      "strandInfo",
      "strand",
      "signatureName",
      "mutationType",
      "contribution",
      "study",
    ],
  },

  {
    sourcePath: "Signature/summary.csv",
    table: "signatureSummary",
    description: "signature summary",
    columns: ["species", "profile", "matrix", "signatureSetName", "count"],
  },

  {
    description: "Backfill study for reference signatures",
    type: "postImport",
    callback: async (connection, logger) => {
      // reference/experimental/published signatures have no study; de novo rows already carry
      // their own study and must not be relabelled, so source excludes them explicitly
      const start = Date.now();
      try {
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
  },

  {
    description: "Refresh materialized views",
    type: "postImport",
    callback: async (connection, logger) => {
      const views = ["signatureOption"];

      for (const view of views) {
        const start = Date.now();
        try {
          logger.info(`Refreshing materialized view "${view}"...`);
          await connection.query(`refresh materialized view "${view}"`);
          const duration = ((Date.now() - start) / 1000).toFixed(2);
          logger.info(
            `Successfully refreshed materialized view "${view}" in ${duration}s`,
          );
        } catch (error) {
          const duration = ((Date.now() - start) / 1000).toFixed(2);
          logger.error(
            `Failed to refresh materialized view "${view}" after ${duration}s: ${error.message}`,
          );
          throw error;
        }
      }
    },
  },
];
