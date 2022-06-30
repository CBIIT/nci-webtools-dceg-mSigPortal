import { fileURLToPath, pathToFileURL } from "url";
import { createRequire } from "module";
import minimist from "minimist";
import { SQSClient } from "@aws-sdk/client-sqs";
import { processMessages } from "./services/queue.js";
import { loadAwsCredentials } from "./services/utils.js";
import { importData, createCustomLogger } from "./startDatabaseImport.js";
import { getSourceProvider } from "./importDatabase.js";

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
  const logger = createCustomLogger("msigportal-queue-worker");
  startQueueWorker(config, schema, sources, sourceProvider, logger);
}

export function startQueueWorker(config, schema, sources, sourceProvider, logger) {
  logger.info("Started msigportal queue worker");

  processMessages({
    sqs: new SQSClient(),
    queueName: config.aws.sqsName,
    visibilityTimeout: config.aws.sqsVisibilityTimeout || 300,
    pollInterval: config.aws.sqsPollInterval || 5,
    messageHandler: async (message) => {
      logger.info("Retrieved message from SQS queue");
      logger.info(message);

      switch (message.type) {
        case "importData":
          return await importData(config, schema, sources, sourceProvider, logger);
        default:
          logger.error(`Unknown message type`);
          return false;
      }
    },
  });
}
