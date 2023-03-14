import { Router } from 'express';
import { validate } from 'uuid';
import path from 'path';
import { downloadDirectory, uploadDirectory } from '../../s3.js';
import logger from '../../logger.js';
import { mkdirs, writeJson, readJson } from '../../utils.js';
import config from '../../../config.json' assert { type: 'json' };

export async function submit(req, res, next) {
  const id = req.params.id;
  if (!validate(id)) res.status(500).json('Invalid ID');

  const inputFolder = path.resolve(config.folders.input, id);
  const outputFolder = path.resolve(config.folders.output, id);
  const paramsFilePath = path.resolve(inputFolder, 'params.json');
  const statusFilePath = path.resolve(outputFolder, 'status.json');
  await mkdirs([inputFolder, outputFolder]);

  const status = {
    id,
    status: 'SUBMITTED',
    submittedAt: new Date(),
  };

  await writeJson(paramsFilePath, req.body);
  await writeJson(statusFilePath, status);

  const s3ClientConfig = { region: config.aws.region };

  await uploadDirectory(
    inputFolder,
    path.join(config.aws.inputKeyPrefix, id),
    config.aws.bucket,
    s3ClientConfig
  );

  await uploadDirectory(
    outputFolder,
    path.join(config.aws.outputKeyPrefix, id),
    config.aws.bucket,
    s3ClientConfig
  );

  try {
    fetch(`${config.email.baseUrl}/extraction/run/${id}`);
  } catch (error) {
    next(error);
  }

  res.json(status);
}

// downloads latest files from s3 and returns status, params, and manifest
async function getJobStatus(id) {
  if (!validate(id)) return `${id} is not a valid ID`;
  const inputFolder = path.resolve(config.folders.input, id);
  const outputFolder = path.resolve(config.folders.output, id);

  try {
    await downloadDirectory(
      inputFolder,
      path.join(config.aws.inputKeyPrefix, id),
      config.aws.bucket,
      { region: config.aws.region }
    );
    await downloadDirectory(
      outputFolder,
      path.join(config.aws.outputKeyPrefix, id),
      config.aws.bucket,
      { region: config.aws.region }
    );

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
  try {
    const id = req.params.id;
    const data = await getJobStatus(id);
    res.json(data);
  } catch (error) {
    logger.error('/refreshExtraction Error');
    next(error);
  }
}

export async function refreshMulti(req, res, next) {
  try {
    const ids = req.body;
    const stauses = await Promise.all(ids.map(getJobStatus));
    res.json(stauses);
  } catch (error) {
    logger.error('/refreshExtractionMulti Error');
    next(error);
  }
}

const router = Router();
router.post('/submitExtraction/:id?', submit);
router.get('/refreshExtraction/:id?', refresh);
router.post('/refreshExtractionMulti', refreshMulti);

export { router };
