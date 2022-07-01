import { fileURLToPath, pathToFileURL } from "url";
import { createRequire } from "module";
import { format } from "util";
import minimist from "minimist";
import { getLogger } from "./services/logger.js";
import { CustomTransport } from "./services/transports.js";
import { loadAwsCredentials, createConnection, createPostgresClient } from "./services/utils.js";
import { importDatabase, getSourceProvider } from "./importDatabase.js";

// determine if this script was launched from the command line
const isMainModule = process.argv[1] === fileURLToPath(import.meta.url);
const require = createRequire(import.meta.url);

if (isMainModule) {
  const config = require("./config.json");
  loadAwsCredentials(config.aws);

  const args = minimist(process.argv.slice(2));
  const schemaPath = pathToFileURL(args.schema || "./schema.js");
  const sourcesPath = pathToFileURL(args.sources || "./sources.js");

  const providerName = args.provider || "s3";
  const defaultProviderArgs = {
    local: ["."],
    s3: [`s3://${config.aws.s3DataBucket}/${config.aws.s3AnalysisKey}`],
  }[providerName];
  const providerArgs = args._.length ? args._ : defaultProviderArgs;

  const { schema } = await import(schemaPath);
  const { sources } = await import(sourcesPath);
  const sourceProvider = getSourceProvider(providerName, providerArgs);
  const logger = createCustomLogger("msigportal-data-import");
  await importData(config, schema, sources, sourceProvider, logger);
  process.exit(0);
}

export async function importData(config, schema, sources, sourceProvider, logger) {
  const connection = createConnection(config.database);
  const logConnection = await createPostgresClient(config.database); // use a separate connection for all log operations
  const importLog = await getPendingImportLog(connection);
  const forceRecreate = importLog.type === "full";

  logger.info(`Started msigportal data import: ${importLog.type}`);

  async function updateImportLog(params) {
    await connection("importLog")
      .where({ id: importLog.id })
      .update({ ...params, updatedAt: new Date() })
      .connection(logConnection);
  }

  async function handleLogEvent(event) {
    const logMessage = `${event.timestamp} ${format(event.message)}`;
    const log = connection.raw(`concat("log", '\n', ?::text)`, [logMessage]);
    await updateImportLog({ log });
  }

  async function shouldCancelImport() {
    const results = await connection("importLog")
      .where({ id: importLog.id, status: "CANCELLED" })
      .connection(logConnection);
    return results.length > 0;
  }

  try {
    logger.customTransport.setHandler(handleLogEvent);
    await updateImportLog({ status: "IN PROGRESS" });
    await importDatabase(connection, schema, sources, sourceProvider, logger, forceRecreate, shouldCancelImport);
    await updateImportLog({ status: "COMPLETED" });
  } catch (exception) {
    logger.error(exception.stack);
    await updateImportLog({ status: "FAILED" });
  } finally {
    logger.customTransport.setHandler(null);
  }

  return true;
}

export function createCustomLogger(name) {
  const logger = getLogger(name);
  logger.customTransport = new CustomTransport();
  logger.add(logger.customTransport);
  return logger;
}

export async function getPendingImportLog(connection) {
  const pendingImportLog = await connection("importLog")
    .where({ status: "PENDING" })
    .orderBy("createdAt", "asc")
    .first();
  return pendingImportLog || (await createImportLog(connection));
}

export async function createImportLog(connection) {
  await connection("importLog").insert({ status: "PENDING" });
  return await getPendingImportLog(connection);
}