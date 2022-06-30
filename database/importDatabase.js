import { basename, resolve } from "path";
import { fileURLToPath, pathToFileURL } from "url";
import { createRequire } from "module";
import minimist from "minimist";
import {
  createConnection,
  createRecordIterator,
  getFileMetadataFromPath,
  importTable,
  initializeSchemaForImport,
  loadAwsCredentials,
  withDuration,
} from "./services/utils.js";
import { LocalProvider } from "./services/providers/localProvider.js";
import { S3Provider } from "./services/providers/s3Provider.js";
import { getLogger } from "./services/logger.js";
import { S3Client } from "@aws-sdk/client-s3";

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
  const forceRecreate = args.recreate;
  loadAwsCredentials(aws);

  const connection = createConnection(database);
  const { schema } = await import(schemaPath);
  const { sources } = await import(sourcesPath);
  const sourceProvider = getSourceProvider(providerName, providerArgs);
  const logger = getLogger("import");
  await importDatabase(connection, schema, sources, sourceProvider, logger, forceRecreate);
  process.exit(0);
}

export async function importDatabase(
  connection,
  schema,
  sources,
  sourceProvider,
  logger,
  forceRecreate = false,
  shouldCancel = async () => false
) {
  const tableSources = sources.filter((source) => source.table);
  const postImportSources = sources.filter((source) => source.type === "postImport");
  let totalCount = 0;

  let { results, duration } = await withDuration(async () => {
    return await connection.transaction(async (connection) => {
      await initializeSchemaForImport(connection, schema, sources, forceRecreate);

      for (let source of tableSources) {
        const { description, table, columns, skipImport, parseConfig } = source;
        const shouldSkip = skipImport ? await skipImport(connection) : () => false;
        const sourcePaths = await getSourcePaths(source, sourceProvider);

        for (const sourcePath of sourcePaths) {
          if (await shouldCancel()) {
            throw new Error(`Cancelled import`);
          }

          const metadata = getFileMetadataFromPath(sourcePath);
          if (shouldSkip(metadata)) continue;

          logger.info(`Importing ${sourcePath} => ${table} (${description})`);
          const { results, duration } = await withDuration(async () => {
            const records = await createRecordIterator(sourcePath, sourceProvider, { columns, parseConfig, logger });
            return await importTable(connection, records, table, logger);
          });

          logger.info(getStatusMessage({ results, duration }));
          totalCount += results;
        }
      }
      return totalCount;
    });
  });

  logger.info(getStatusMessage({ results, duration }));
}

function getStatusMessage({ results, duration }) {
  return `Finished importing ${results} rows in ${duration.toFixed(2)}s (${Math.round(results / duration)} rows/s)`;
}

async function getSourcePaths(source, sourceProvider) {
  return source.type === "folder" ? await sourceProvider.listFiles(source.sourcePath) : [source.sourcePath];
}

export function getSourceProvider(providerName, providerArgs) {
  switch (providerName) {
    case "local":
      return new LocalProvider(...providerArgs);
    case "s3":
      return new S3Provider(new S3Client(), ...providerArgs);
    default:
      throw new Error(`Unknown provider: ${providerName}`);
  }
}
