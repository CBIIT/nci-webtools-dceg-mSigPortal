import path from 'path';
import { readdir, stat, unlink, rmdir } from 'fs/promises';
import { CronJob } from 'cron';

const FOURTEEN_DAYS_MS = 14 * 24 * 60 * 60 * 1000;

/**
 * Recursively deletes files and empty directories whose mtime exceeds maxAgeMs.
 * Mirrors `find <dir> -mtime +14 -delete` (depth-first, bottom-up).
 * The root directory itself is never deleted.
 */
async function deleteOldEntries(dir, maxAgeMs, logger, isRoot = true) {
  let entries;
  try {
    entries = await readdir(dir, { withFileTypes: true });
  } catch (err) {
    logger.error(`Failed to read directory ${dir}: ${err.message}`);
    return;
  }

  for (const entry of entries) {
    const fullPath = path.join(dir, entry.name);

    if (entry.isDirectory()) {
      await deleteOldEntries(fullPath, maxAgeMs, logger, false);

      if (isRoot) continue;

      try {
        const dirStat = await stat(fullPath);
        if (Date.now() - dirStat.mtimeMs > maxAgeMs) {
          await rmdir(fullPath);
          logger.info(`Removed directory: ${fullPath}`);
        }
      } catch (err) {
        if (err.code !== 'ENOTEMPTY' && err.code !== 'ENOENT') {
          logger.error(`Failed to remove directory ${fullPath}: ${err.message}`);
        }
      }
    } else {
      try {
        const fileStat = await stat(fullPath);
        if (Date.now() - fileStat.mtimeMs > maxAgeMs) {
          await unlink(fullPath);
          logger.info(`Removed file: ${fullPath}`);
        }
      } catch (err) {
        if (err.code !== 'ENOENT') {
          logger.error(`Failed to remove file ${fullPath}: ${err.message}`);
        }
      }
    }
  }
}

async function cleanup(logger) {
  const { INPUT_FOLDER, OUTPUT_FOLDER } = process.env;
  await deleteOldEntries(path.resolve(INPUT_FOLDER), FOURTEEN_DAYS_MS, logger);
  await deleteOldEntries(path.resolve(OUTPUT_FOLDER), FOURTEEN_DAYS_MS, logger);
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
