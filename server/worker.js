import path from 'path';
import knex from 'knex';
import { isMainModule, readJson } from './services/utils.js';
import { createLogger } from './services/logger.js';
import { mkdirs } from './services/utils.js';
import { profilerExtraction } from './services/api/visualization/profilerExtraction.js';

if (isMainModule(import.meta)) {
  try {
    await main(process.argv, process.env);
    process.exit(0);
  } catch (e) {
    console.error(e);
    process.exit(1);
  }
}

export async function main(argv = process.argv, env = process.env) {
  const id = argv[2];
  if (!id) throw new Error('Missing id');

  const inputFolder = path.resolve(env.INPUT_FOLDER, id);
  const outputFolder = path.resolve(env.OUTPUT_FOLDER, id);
  await mkdirs([inputFolder, outputFolder]);

  const paramsFilePath = path.resolve(inputFolder, 'params.json');
  const params = await readJson(paramsFilePath);
  const logger = createLogger(env.APP_NAME, env.LOG_LEVEL);
  const localDb = knex({
    client: 'better-sqlite3',
    connection: {
      filename: path.join(outputFolder, `local.sqlite3`),
    },
    useNullAsDefault: true,
  });
  logger.log({ params });
  return await profilerExtraction(params, logger, localDb, env);
}
