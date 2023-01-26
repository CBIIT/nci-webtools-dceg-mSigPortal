import path from 'path';
import { mkdirs, writeJson } from './utils.js';
import { getSqliteConnection } from './database.js';
import { getWorker } from './workers.js';
import { existsSync } from 'fs';
const { WORKER_TYPE } = process.env;

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
