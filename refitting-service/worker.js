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
  
  // Log what we read from params.json
  logger.info(`[${id}] Reading params.json from: ${paramsFilePath}`);
  logger.info(`[${id}] params.json contents:`, JSON.stringify(params, null, 2));
  logger.info(`[${id}] params.email = ${params.email}`);
  logger.info(`[${id}] params.jobName = ${params.jobName}`);

  // Get file paths from input directory
  const fs = await import('fs');
  const files = fs.readdirSync(inputFolder);
  const mafFile = files.find(f => f.startsWith('mafFile_'));
  const genomicFile = files.find(f => f.startsWith('genomicFile_'));
  const clinicalFile = files.find(f => f.startsWith('clinicalFile_'));

  if (!mafFile || !genomicFile || !clinicalFile) {
    throw new Error('Missing required input files (mafFile, genomicFile, or clinicalFile)');
  }

  // Construct complete params object with file paths
  const completeParams = {
    ...params,
    id,
    mafFile: path.resolve(inputFolder, mafFile),
    genomicFile: path.resolve(inputFolder, genomicFile),
    clinicalFile: path.resolve(inputFolder, clinicalFile),
    outputPath: outputFolder,
  };

  logger.info(`[${id}] Complete params being passed to refitting:`, JSON.stringify(completeParams, null, 2));
  logger.info(`[${id}] completeParams.email = ${completeParams.email}`);
  logger.info(`[${id}] completeParams.jobName = ${completeParams.jobName}`);
  
  logger.log({ params: completeParams });
  
  return await refitting(completeParams, logger, env);
}
