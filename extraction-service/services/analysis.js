import path from 'path';
import { mkdirs, writeJson } from './utils.js';
import { getSqliteConnection } from './database.js';
import { getWorker } from './workers.js';
import { existsSync } from 'fs';
const { WORKER_TYPE } = process.env;
import { extraction } from './extraction.js';
import { downloadDirectory } from './s3.js';
import { readJson } from './utils.js';

export async function submit(params, app, env = process.env) {
  const id = params.id;
  const inputFolder = path.resolve(env.INPUT_FOLDER, id);
  const outputFolder = path.resolve(env.OUTPUT_FOLDER, id);
  const paramsFilePath = path.resolve(inputFolder, 'params.json');
  const statusFilePath = path.resolve(outputFolder, 'status.json');
  await mkdirs([inputFolder, outputFolder]);

  const worker = getWorker(WORKER_TYPE);
  const status = { id, status: 'SUBMITTED' };

  await writeJson(paramsFilePath, params);
  await writeJson(statusFilePath, status);

  worker(id, app).catch(console.error);
  return status;
}

export async function query(params, env = process.env) {
  const { id, table, columns, conditions, orderBy, offset, limit } = params;
  const databaseFilePath = path.resolve(env.OUTPUT_FOLDER, id, 'results.db');
  if (!existsSync(databaseFilePath)) return [];
  return await getSqliteConnection(databaseFilePath)
    .select(columns || '*')
    .from(table)
    .where(conditions || {})
    .offset(offset || 0)
    .limit(limit || 100000)
    .orderBy(orderBy || []);
}

export async function run(id, app, env = process.env) {
  if (!id) throw new Error('Missing id');
  app.locals.logger.debug(id);
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
  const logger = app.locals.logger;
  const dbConnection = app.locals.connection;
  logger.log({ params });
  return await extraction(params, logger, dbConnection, env);
}
