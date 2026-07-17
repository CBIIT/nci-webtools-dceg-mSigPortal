import { fileURLToPath, pathToFileURL } from "url";
import minimist from "minimist";
import { createLogger, formatObject } from "./services/logger.js";
import {
  getConfigFromEnv,
  loadAwsCredentials,
  createConnection,
  getSourceProvider,
  registerKnexConnection,
  clearKnexConnection,
  abortActiveConnections,
} from "./services/utils.js";
import { sendImportNotification } from "./services/notifications.js";
import { importDatabase } from "./importDatabase.js";

// determine if this script was launched from the command line
const isMainModule = process.argv[1] === fileURLToPath(import.meta.url);

let shuttingDown = false;

if (isMainModule) {
  const config = getConfigFromEnv();
  loadAwsCredentials(config.aws);

  const args = minimist(process.argv.slice(2));
  const schemaPath = pathToFileURL(args.schema || "./schema.js");
  const sourcesPath = pathToFileURL(args.sources || "./sources.js");

  const providerName = args.provider || "s3";
  const defaultProviderArgs = {
    local: ["."],
    s3: [`s3://${config.data.bucket}/${config.data.s3}`],
  }[providerName];
  const providerArgs = args._.length ? args._ : defaultProviderArgs;

  const { schema } = await import(schemaPath);
  const { sources } = await import(sourcesPath);
  const sourceProvider = getSourceProvider(providerName, providerArgs);
  const logger = createLogger(
    process.env.APP_NAME ? `${process.env.APP_NAME}-data-import` : "msigportal-data-import",
    process.env.LOG_LEVEL || "info"
  );

  const importPromise = importData(
    config,
    schema,
    sources,
    sourceProvider,
    logger
  );

  async function handleShutdown(signal) {
    if (shuttingDown) {
      return;
    }
    shuttingDown = true;
    logger.warn(`Received ${signal}; aborting import and sending notification`);
    await abortActiveConnections();
    try {
      await importPromise;
    } catch {
      // importData marks FAILED and sends email in finally
    }
    process.exit(1);
  }

  process.on("SIGTERM", () => {
    handleShutdown("SIGTERM").catch((error) => {
      logger.error(`Shutdown handler failed: ${error.stack}`);
      process.exit(1);
    });
  });
  process.on("SIGINT", () => {
    handleShutdown("SIGINT").catch((error) => {
      logger.error(`Shutdown handler failed: ${error.stack}`);
      process.exit(1);
    });
  });

  try {
    await importPromise;
    process.exit(0);
  } catch (exception) {
    logger.error(exception.stack);
    process.exit(1);
  }
}

/**
 * Wraps a winston logger so each message is also collected for the import email body.
 */
export function createCapturingLogger(logger, logLines) {
  const capture =
    (level) =>
    (message, ...rest) => {
      const timestamp = new Date().toISOString();
      logLines.push(`${timestamp} [${level}] ${formatObject(message)}`);
      return logger[level](message, ...rest);
    };

  return {
    info: capture("info"),
    warn: capture("warn"),
    error: capture("error"),
    debug: capture("debug"),
  };
}

export async function importData(
  config,
  schema,
  sources,
  sourceProvider,
  logger
) {
  const connection = createConnection(config.database);
  registerKnexConnection(connection);
  const importLog = await getPendingImportLog(connection);
  let status = "COMPLETED";
  const logLines = [];
  const capturingLogger = createCapturingLogger(logger, logLines);

  capturingLogger.info(`Started msigportal data import`);

  async function updateImportLog(params) {
    await connection("importLog")
      .where({ id: importLog.id })
      .update({ ...params, updatedAt: new Date() });
  }

  async function shouldCancelImport() {
    if (shuttingDown) {
      return true;
    }
    const results = await connection("importLog").where({
      id: importLog.id,
      status: "CANCELLED",
    });
    return results.length > 0;
  }

  try {
    await updateImportLog({ status: "IN PROGRESS" });
    await importDatabase(
      config.database,
      schema,
      sources,
      sourceProvider,
      capturingLogger,
      shouldCancelImport
    );
    if (shuttingDown) {
      throw new Error("Import cancelled due to shutdown signal");
    }
    await updateImportLog({ status: "COMPLETED", log: logLines.join("\n") });
  } catch (exception) {
    status = "FAILED";
    capturingLogger.error(exception.stack);
    try {
      await updateImportLog({ status: "FAILED", log: logLines.join("\n") });
    } catch (updateError) {
      capturingLogger.error(
        `Failed to update import log status: ${updateError.message}`
      );
    }
    throw exception;
  } finally {
    try {
      await sendImportNotification({
        status,
        startTime: new Date(importLog.createdAt).toString(),
        logs: logLines.join("\n") || null,
        env: process.env,
      });
    } catch (exception) {
      logger.error(`Failed to send import notification: ${exception.stack}`);
    } finally {
      clearKnexConnection(connection);
      try {
        await connection.destroy();
      } catch {
        // ignore destroy errors
      }
    }
  }

  return true;
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
