import path from 'path';
import { readdir, stat, rm } from 'fs/promises';
import { CronJob } from 'cron';

const MAX_AGE_MS = 14 * 24 * 60 * 60 * 1000;

/**
 * Removes immediate child entries of dir whose mtime exceeds MAX_AGE_MS.
 * Directories are removed recursively (rm -rf).
 */
async function removeOldEntries(dir, logger) {
  let entries;
  try {
    entries = await readdir(dir, { withFileTypes: true });
  } catch (err) {
    logger.error(`Failed to read directory ${dir}: ${err.message}`);
    return;
  }

  for (const entry of entries) {
    const fullPath = path.join(dir, entry.name);
    try {
      const { mtimeMs } = await stat(fullPath);
      if (Date.now() - mtimeMs > MAX_AGE_MS) {
        await rm(fullPath, { recursive: true, force: true });
        logger.info(`Removed: ${fullPath}`);
      }
    } catch (err) {
      if (err.code !== 'ENOENT') {
        logger.error(`Failed to remove ${fullPath}: ${err.message}`);
      }
    }
  }
}

async function cleanup(logger) {
  const { INPUT_FOLDER, OUTPUT_FOLDER } = process.env;
  await removeOldEntries(path.resolve(INPUT_FOLDER), logger);
  await removeOldEntries(path.resolve(OUTPUT_FOLDER), logger);
}

export function startCron(app) {
  const { logger } = app.locals;
  const job = new CronJob(
    '0 12 * * *',
    async () => {
      logger.info('Starting cleanup of expired jobs');
      try {
        await cleanup(logger);
        logger.info('Cleanup complete');
      } catch (err) {
        logger.error(`Cleanup failed: ${err.message}`);
      }
    },
    null,
    true,
    'America/New_York'
  );
  logger.info('Cron job scheduled: cleanup expired data daily at 12:00 PM ET');
  return job;
}
