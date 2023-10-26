import { fileURLToPath, pathToFileURL } from "url";
import { createRequire } from "module";
import minimist from "minimist";
import {
  createPostgresConnection,
  initializeSchemaForImport,
  importPostgresTable,
  loadAwsCredentials,
  withDuration,
  getSourceProvider,
} from "./services/utils.js";
import { getLogger } from "./services/logger.js";

// determine if this script was launched from the command line
const isMainModule = process.argv[1] === fileURLToPath(import.meta.url);
const require = createRequire(import.meta.url);

if (isMainModule) {
  const { aws, database } = require("./config.json");
  const args = minimist(process.argv.slice(2));
  const schemaPath = pathToFileURL(args.schema || "./schema.js");
  const sourcesPath = pathToFileURL(args.sources || "./sources.js");
  const providerName = args.provider || "local";
  const providerArgs = [...args._];
  loadAwsCredentials(aws);

  const { schema } = await import(schemaPath);
  const { sources } = await import(sourcesPath);
  const sourceProvider = getSourceProvider(providerName, providerArgs);
  const logger = getLogger("import");
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
  const postImportSteps = sources.filter((source) => source.type === "postImport");
  const connection = await createPostgresConnection(connectionConfig);
  let totalCount = 0;

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
        return await importPostgresTable(connection, inputStream, table, columns);
      });

      totalCount += results;
      logger.info(getStatusMessage({results, duration}));      
    }

    for (let postImportStep of postImportSteps) {
      logger.info(`Running post-import step (${postImportStep.description})`);
      await postImportStep.callback(connection);
    }

    return totalCount;
  });

  logger.info(getStatusMessage({results, duration}));      
}

function getStatusMessage({ results, duration }) {
  return `Finished importing ${results} rows in ${duration.toFixed(2)}s (${Math.round(results / duration)} rows/s)`;
}
