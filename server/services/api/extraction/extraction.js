import { Router } from 'express';
import { validate } from 'uuid';
import path from 'path';
import { mkdirs, writeJson, readJson } from '../../utils.js';
import { getWorker } from '../../workers.js';
import { exampleProcessor } from './exampleProcessor.js';
import fs from 'fs';
import { copy } from 'fs-extra';
import { randomUUID } from 'crypto';
import { handleValidationErrors } from '../../middleware.js';
import { check } from 'express-validator';
const env = process.env;

export async function submit(req, res, next) {
  const id = req.params.id;
  if (!validate(id)) res.status(500).json('Invalid ID');

  const inputFolder = path.resolve(env.INPUT_FOLDER, id);
  const outputFolder = path.resolve(env.OUTPUT_FOLDER, id);
  const paramsFilePath = path.resolve(inputFolder, 'params.json');
  const statusFilePath = path.resolve(outputFolder, 'status.json');
  await mkdirs([inputFolder, outputFolder]);

  const status = {
    id,
    status: 'SUBMITTED',
    submittedAt: new Date(),
  };

  const worker = getWorker(env.EXTRACTION_WORKER_TYPE || 'local');

  await writeJson(paramsFilePath, req.body);
  await writeJson(statusFilePath, status);

  worker(id, req.app, 'extraction', env);
  res.json(status);
}

// downloads latest files from s3 and returns status, params, and manifest
async function getJobStatus(id) {
  //if (!validate(id)) return `${id} is not a valid ID`;
  if (!/^Example_/.test(id) && !validate(id)) {
    return `${id} is not a valid ID`;
  }
  try {
    const inputFolder = path.resolve(env.INPUT_FOLDER, id);
    const outputFolder = path.resolve(env.OUTPUT_FOLDER, id);

    const paramsFilePath = path.resolve(inputFolder, 'params.json');
    const statusFilePath = path.resolve(outputFolder, 'status.json');
    const manifestFilePath = path.resolve(outputFolder, 'manifest.json');
    const params = await readJson(paramsFilePath);
    const status = await readJson(statusFilePath);
    const manifest = await readJson(manifestFilePath);

    return { params, status, manifest };
  } catch (error) {
    return error.message;
  }
}

export async function refresh(req, res, next) {
  const { logger } = req.app.locals;
  try {
    const id = req.params.id;
    const data = await getJobStatus(id);
    if (typeof data === 'string') {
      res.status(500).json(data);
    } else {
      res.json(data);
    }
  } catch (error) {
    logger.error('/refreshExtraction Error');
    next(error);
  }
}

export async function refreshMulti(req, res, next) {
  const { logger } = req.app.locals;
  try {
    const ids = req.body;
    const stauses = await Promise.all(ids.map(getJobStatus));
    res.json(stauses);
  } catch (error) {
    logger.error('/refreshExtractionMulti Error');
    next(error);
  }
}

export async function extractionExample(req, res, next) {
  const { logger } = req.app.locals;
  try {
    const id = req.params.id;
    const uuid = randomUUID();
    const dataFolder = path.resolve(env.DATA_FOLDER);
    const outputFolder = path.resolve(env.OUTPUT_FOLDER, uuid);
    const exampleFolderPath = path.resolve(
      dataFolder,
      'examples',
      'extraction',
      id
    );

    if (fs.existsSync(exampleFolderPath)) {
      const exampleOutputFolderName = path.basename(outputFolder);
      //copy example folder into outputFolder
      await copy(exampleFolderPath, outputFolder);
      const result = await exampleProcessor(
        exampleOutputFolderName,
        id,
        uuid,
        env
      );
      res.json(result);
    } else {
      res.status(404).json({
        error: 'Example folder does not exist: ',
        exampleFolderPath,
      });
    }
  } catch (error) {
    logger.error('/extractionExample Error', error);
    next(error);
  }
}

// validate parameter limits for extraction submit
const minMaxValidation = (param, min, max) =>
  check(param)
    .custom((e) => e >= min && e <= max)
    .withMessage(`${param} must be between ${min} and ${max}`);

const router = Router();
router.post(
  '/submitExtraction/:id?',
  minMaxValidation('args.minimum_signatures', 1, 20),
  minMaxValidation('args.maximum_signatures', 1, 20),
  minMaxValidation('args.nmf_replicates', 1, 200),
  minMaxValidation('args.min_nmf_iterations', 1, 20000),
  minMaxValidation('args.max_nmf_iterations', 1, 2000000),
  minMaxValidation('args.nmf_test_conv', 1, 20000),
  handleValidationErrors,
  submit
);
router.get('/refreshExtraction/:id?', refresh);
router.post('/refreshExtractionMulti', refreshMulti);
router.get('/extractionExample/:id', extractionExample);

export { router };
