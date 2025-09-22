import { createLogger } from "./logger.js";
import rWrapper from "r-wrapper";
const r = rWrapper.async;

const [jobId] = process.argv.slice(2);
const logger = createLogger(`Job ${jobId}`, process.env?.LOG_LEVEL || "debug");

if (!jobId) {
  logger.error("No jobId specified");
  process.exit(1);
}

try {
  await refitting(jobId);
  process.exit(0);
} catch (e) {
  logger.error(e);
  process.exit(1);
}

async function refitting(jobId) {
  const start = new Date();
  try {
    const sum = await r("./refitting.R", "sum", { a: 1, b: 2 }); // test call
    logger.info(`R test call result: ${sum}`);
    const end = new Date();
    const duration = (end - start) / 1000;
    logger.info(`Refitting job ${jobId} completed in ${duration} seconds`);
  } catch (e) {
    const done = new Date();
    const duration = (done - start) / 1000;
    logger.info(`Refitting job ${jobId} failed after ${duration} seconds`);
    logger.error(e);
    throw e;
  }
}
