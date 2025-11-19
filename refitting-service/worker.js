import path from "path";
import { isMainModule, readJson, mkdirs } from "./services/utils.js";
import { createLogger } from "./services/logger.js";
import { refitting } from "./services/refitting.js";
import fs from "fs";

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
  const logger = createLogger(`Refitting Worker ${id}`, env.LOG_LEVEL);
  if (!id) throw new Error("Missing id");

  const inputFolder = path.resolve(env.INPUT_FOLDER, id);
  const outputFolder = path.resolve(env.OUTPUT_FOLDER, id);
  await mkdirs([inputFolder, outputFolder]);

  const paramsFilePath = path.resolve(inputFolder, "params.json");
  const params = await readJson(paramsFilePath);

  const mafFile = path.resolve(inputFolder, params.mafFileName);
  const genomicFile = path.resolve(inputFolder, params.genomicFileName);
  const clinicalFile = path.resolve(inputFolder, params.clinicalFileName);

  // Validate files exist
  if (!fs.existsSync(mafFile)) {
    logger.error(`MAF file not found: ${mafFile}`);
    throw new Error(`MAF file not found: ${mafFile}`);
  }
  if (!fs.existsSync(genomicFile)) {
    logger.error(`Genomic file not found: ${genomicFile}`);
    throw new Error(`Genomic file not found: ${genomicFile}`);
  }
  if (!fs.existsSync(clinicalFile)) {
    logger.error(`Clinical file not found: ${clinicalFile}`);
    throw new Error(`Clinical file not found: ${clinicalFile}`);
  }

  // Construct complete params object with file paths
  const completeParams = {
    ...params,
    id,
    mafFile,
    genomicFile,
    clinicalFile,
    outputPath: outputFolder,
    profileType: params.signatureType,
  };

  logger.info({ params: completeParams });
  return await refitting(completeParams, logger, env);
}
