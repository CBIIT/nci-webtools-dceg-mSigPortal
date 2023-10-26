import { Router } from 'express';
import { validate } from 'uuid';
import path from 'path';
import { parseCSV } from '../general.js';
import { schema } from './userSchema.js';
import { getSignatureData } from '../../query.js';
import { mkdirs } from '../../utils.js';
import { sqliteImport } from '../../sqlite.js';

const env = process.env;

async function submit(req, res, next) {
  const { logger } = req.app.locals;
  const id = req.params.id;
  if (!validate(id)) res.status(500).json('Invalid ID');
  const inputFolder = path.resolve(env.INPUT_FOLDER, id);
  const outputFolder = path.resolve(env.OUTPUT_FOLDER, id);
  await mkdirs([inputFolder, outputFolder]);
  const { exposureFile, matrixFile, signatureFile, signatureSetName } =
    req.body;
  const exposurePath = path.resolve(inputFolder, exposureFile);
  const matrixPath = path.resolve(inputFolder, matrixFile);
  const signaturePath = signatureFile
    ? path.resolve(inputFolder, signatureFile)
    : '';
  const exposureData = await parseCSV(exposurePath);
  const matrixData = await parseCSV(matrixPath);
  const signatureData = signaturePath ? await parseCSV(signaturePath) : '';

  // transform input data into format suitable for db import
  const transformExposure = exposureData
    .map((e) => {
      const keys = Object.keys(e);
      const sampleName = keys[0];
      const signatureNames = keys.slice(1);
      return signatureNames.map((signatureName) => ({
        signatureName,
        sample: e[sampleName],
        exposure: e[signatureName],
        cancer: 'Input',
      }));
    })
    .flat();
  const transformMatrix = matrixData
    .map((e) => {
      const keys = Object.keys(e);
      const mutationType = keys[0];
      const sampleNames = keys.slice(1);
      return sampleNames.map((sample) => ({
        sample,
        mutationType: e[mutationType],
        mutations: e[sample],
        cancer: 'Input',
      }));
    })
    .flat();
  const transformSignature = signatureData
    ? signatureData
        .map((e) => {
          const keys = Object.keys(e);
          const mutationType = keys[0];
          const signatureNames = keys.slice(1);
          return signatureNames.map((signatureName) => ({
            signatureName,
            mutationType: e[mutationType],
            contribution: e[signatureName],
          }));
        })
        .flat()
    : await getSignatureData(
        req.app.locals.connection,
        { signatureSetName },
        ['signatureName', 'mutationType', 'contribution'],
        false
      );

  // import data into user session table
  const connection = req.app.locals.sqlite(id, 'local');
  try {
    const importStatus = await sqliteImport(
      connection,
      {
        exposure: transformExposure,
        seqmatrix: transformMatrix,
        ...(transformSignature && { signature: transformSignature }),
      },
      schema
    );
    if (!importStatus)
      res.status(500).json('Failed to import data into database');
    res.json(id);
  } catch (error) {
    logger.error(error);
    next(error);
  }
}

const router = Router();
router.post('/submitExploration/:id', submit);

export { router };
