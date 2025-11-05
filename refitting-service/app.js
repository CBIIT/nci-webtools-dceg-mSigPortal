import { createLogger } from "./logger.js";
import { refitting } from './services/refitting.js';
import path from 'path';
import fs from 'fs';

// Parse command line arguments
const args = process.argv.slice(2);
const jobId = args[0];

if (!jobId) {
  console.error("No jobId specified");
  process.exit(1);
}

// Parse additional arguments
const argMap = {};
for (let i = 1; i < args.length; i += 2) {
  if (args[i].startsWith('--') && i + 1 < args.length) {
    argMap[args[i].substring(2)] = args[i + 1];
  }
}

const logger = createLogger(`Job ${jobId}`, process.env?.LOG_LEVEL || "debug");

try {
  // Add jobId to params
  const params = {
    ...argMap,
    id: jobId
  };
  
  await refitting(params, logger, process.env);
  process.exit(0);
} catch (e) {
  logger.error(e);
  process.exit(1);
}