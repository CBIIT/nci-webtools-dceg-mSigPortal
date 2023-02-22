import path from 'path';
import knex from 'knex';
import { isMainModule, readJson } from './services/utils.js';
import { createLogger } from './services/logger.js';
import { extraction } from './services/extraction.js';
import { downloadDirectory } from './services/s3.js';
import { mkdirs } from './services/utils.js';

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

  // download folders from s3
  await downloadDirectory(
    inputFolder,
    path.join(env.INPUT_KEY_PREFIX, id),
    env.DATA_BUCKET,
    { region: env.AWS_DEFAULT_REGION }
  );
  await downloadDirectory(
    outputFolder,
    path.join(env.OUTPUT_KEY_PREFIX, id),
    env.DATA_BUCKET,
    { region: env.AWS_DEFAULT_REGION }
  );

  const paramsFilePath = path.resolve(inputFolder, 'params.json');
  const params = await readJson(paramsFilePath);
  const logger = createLogger(env.APP_NAME, env.LOG_LEVEL);
  const dbConnection = knex({
    client: 'postgres',
    connection: {
      host: env.POSTGRES_HOST,
      port: env.POSTGRES_PORT,
      user: env.POSTGRES_USER,
      password: env.POSTGRES_PASS,
      database: env.POSTGRES_DB,
    },
  });
  logger.log({ params });
  return await extraction(params, logger, dbConnection, env);
}
