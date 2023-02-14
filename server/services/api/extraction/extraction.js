import { Router } from 'express';
import { validate } from 'uuid';
import path from 'path';
import AWS from 'aws-sdk';
import tar from 'tar';
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

  await new AWS.S3()
    .upload({
      Body: tar
        .c({ sync: true, gzip: true, C: config.folders.input }, [id])
        .read(),
      Bucket: config.aws.bucket,
      Key: `${config.aws.inputKeyPrefix}${id}/${id}.tgz`,
    })
    .promise();

  await new AWS.S3()
    .upload({
      Body: tar
        .c({ sync: true, gzip: true, C: config.folders.output }, [id])
        .read(),
      Bucket: config.aws.bucket,
      Key: `${config.aws.outputKeyPrefix}${id}/${id}.tgz`,
    })
    .promise();

  res.json(status);
}

const router = Router();
router.post('/submitExtraction/:id', submit);

export { router };
