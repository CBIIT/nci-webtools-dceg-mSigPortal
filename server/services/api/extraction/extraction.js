import { Router } from 'express';
import { validate } from 'uuid';
import path from 'path';
import { getDirectory, putDirectory } from '../../s3.js';
import logger from '../../logger.js';
import { mkdirs, writeJson } from '../../utils.js';
import config from '../../../config.json' assert { type: 'json' };

export async function submit(req, res, next) {
  const id = req.params.id;
  if (!validate(id)) res.status(500).json('Invalid ID');

  const inputFolder = path.resolve(config.folders.input, id);
  const outputFolder = path.resolve(config.folders.output, id);
  const paramsFilePath = path.resolve(inputFolder, 'params.json');
  const statusFilePath = path.resolve(outputFolder, 'status.json');
  await mkdirs([inputFolder, outputFolder]);

  const status = { id, status: 'SUBMITTED' };

  await writeJson(paramsFilePath, req.body);
  await writeJson(statusFilePath, status);

  const s3ClientConfig = { region: config.aws.region };

  await putDirectory(
    inputFolder,
    path.join(config.aws.inputKeyPrefix, id),
    config.aws.bucket,
    s3ClientConfig
  );

  await putDirectory(
    outputFolder,
    path.join(config.aws.outputKeyPrefix, id),
    config.aws.bucket,
    s3ClientConfig
  );

  res.json(status);
}

export async function refresh(req, res, next) {
  const id = req.params.id;
  if (!validate(id)) res.status(500).json('Invalid ID');

  const inputFolder = path.resolve(config.folders.input, id);
  const outputFolder = path.resolve(config.folders.output, id);
  console.log(
    inputFolder,
    path.join(config.aws.inputKeyPrefix, id),
    config.aws.bucket,
    { region: config.aws.region }
  );
  try {
    await getDirectory(
      inputFolder,
      path.join(config.aws.inputKeyPrefix, id),
      config.aws.bucket,
      { region: config.aws.region }
    );
    await getDirectory(
      outputFolder,
      path.join(config.aws.outputKeyPrefix, id),
      config.aws.bucket,
      { region: config.aws.region }
    );

    res.json(id);
  } catch (error) {
    logger.error('/refreshExtraction Error')
    next(error);
  }
}

const router = Router();
router.post('/submitExtraction/:id?', submit);
router.get('/refreshExtraction/:id?', refresh);

export { router };
