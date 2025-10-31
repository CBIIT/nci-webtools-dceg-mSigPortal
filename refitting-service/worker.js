import path from 'path';
import { isMainModule, readJson, mkdirs } from './services/utils.js';
import { createLogger } from './services/logger.js';
import { refitting } from './services/refitting.js';

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
  
  if (!params) {
    throw new Error(`Failed to read params.json from ${paramsFilePath}`);
  }

  const logger = createLogger(env.REFITTING_APP_NAME || 'RefittingService', env.LOG_LEVEL);
  
  logger.log({ params });
  
  return await refitting(params, logger, env);
}
