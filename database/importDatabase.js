import { fileURLToPath, pathToFileURL } from "url";
import minimist from "minimist";
import {
  getConfigFromEnv,
  createPostgresConnection,
  initializeSchemaForImport,
  importPostgresTable,
  loadAwsCredentials,
  withDuration,
  getSourceProvider,
  registerPgConnection,
  clearPgConnection,
} from "./services/utils.js";
import { createLogger } from "./services/logger.js";

// determine if this script was launched from the command line
const isMainModule = process.argv[1] === fileURLToPath(import.meta.url);

if (isMainModule) {
  const { aws, database } = getConfigFromEnv();
  const args = minimist(process.argv.slice(2));
  const schemaPath = pathToFileURL(args.schema || "./schema.js");
  const sourcesPath = pathToFileURL(args.sources || "./sources.js");
  const providerName = args.provider || "local";
  const providerArgs = [...args._];
  loadAwsCredentials(aws);

  const { schema } = await import(schemaPath);
  const { sources } = await import(sourcesPath);
  const sourceProvider = getSourceProvider(providerName, providerArgs);
  const logger = createLogger("import");
  await importDatabase(database, schema, sources, sourceProvider, logger);
  process.exit(0);
}

export async function importDatabase(
  connectionConfig,
  schema,
  sources,
  sourceProvider,
  logger,
  shouldCancel = async () => false
) {
  const tableSources = sources.filter((source) => source.table);
  const postImportSteps = sources.filter(
    (source) => source.type === "postImport"
  );
  const connection = await createPostgresConnection(connectionConfig);
  registerPgConnection(connection);
  let totalCount = 0;

  try {
    const { results, duration } = await withDuration(async () => {
      await initializeSchemaForImport(connectionConfig, schema, sources);

      for (let source of tableSources) {
        const { description, table, columns, sourcePath } = source;

        if (await shouldCancel()) {
          throw new Error(`Cancelled import`);
        }

        logger.info(`Importing ${sourcePath} => ${table} (${description})`);

        const { results, duration } = await withDuration(async () => {
          const inputStream = await sourceProvider.readFile(sourcePath);
          return await importPostgresTable(
            connection,
            inputStream,
            table,
            columns
          );
        });

        totalCount += results;
        logger.info(getStatusMessage({ results, duration }));
      }

      for (let postImportStep of postImportSteps) {
        logger.info(`Running post-import step (${postImportStep.description})`);
        await postImportStep.callback(connection, logger);
      }

      return totalCount;
    });

    logger.info(getStatusMessage({ results, duration }));
  } finally {
    clearPgConnection(connection);
    try {
      await connection.end();
    } catch (error) {
      // Expected on SIGTERM when abortActiveConnections already ended the client
      logger.warn(`Postgres connection already closed: ${error.message}`);
    }
  }
}

function getStatusMessage({ results, duration }) {
  return `Finished importing ${results} rows in ${duration.toFixed(
    2
  )}s (${Math.round(results / duration)} rows/s)`;
}
