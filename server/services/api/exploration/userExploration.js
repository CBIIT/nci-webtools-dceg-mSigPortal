import { Router } from 'express';
import { validate } from 'uuid';
import { parseCSV } from '../general.js';
import { schema } from './userSchema.js';
import { getSignatureData } from '../../query.js';
import { mkdirs, resolveWithin, sanitizeFilename } from '../../utils.js';
import { sqliteImport } from '../../sqlite.js';

const env = process.env;

async function submit(req, res, next) {
  const { logger } = req.app.locals;
  const id = req.params.id;
  if (!validate(id)) return res.status(500).json('Invalid ID');
  const inputFolder = resolveWithin(env.INPUT_FOLDER, id);
  const outputFolder = resolveWithin(env.OUTPUT_FOLDER, id);
  await mkdirs([inputFolder, outputFolder]);
  const { exposureFile, matrixFile, signatureFile, signatureSetName } =
    req.body;
  const exposurePath = resolveWithin(
    inputFolder,
    sanitizeFilename(exposureFile)
  );
  const matrixPath = resolveWithin(inputFolder, sanitizeFilename(matrixFile));
  const signaturePath = signatureFile
    ? resolveWithin(inputFolder, sanitizeFilename(signatureFile))
    : '';
  const exposureData = await parseCSV(exposurePath);
  const matrixData = await parseCSV(matrixPath);
  const signatureData = signaturePath ? await parseCSV(signaturePath) : '';

  // validate first-column headers (case-insensitive) before importing
  const normalize = (s) =>
    (s || '')
      .replace(/^\uFEFF/, '')
      .trim()
      .toLowerCase();
  const exposureCols = exposureData.length ? Object.keys(exposureData[0]) : [];
  const matrixCols = matrixData.length ? Object.keys(matrixData[0]) : [];
  if (normalize(exposureCols[0]) !== 'samples')
    return res.status(400).json({
      error: `The exposure/activity file's first column must be "Samples" (found "${
        exposureCols[0] ?? ''
      }").`,
    });
  if (exposureCols.length < 2)
    return res.status(400).json({
      error:
        'The exposure/activity file must include at least one signature column.',
    });
  if (normalize(matrixCols[0]) !== 'mutationtype')
    return res.status(400).json({
      error: `The mutation matrix file's first column must be "MutationType" (found "${
        matrixCols[0] ?? ''
      }").`,
    });
  if (matrixCols.length < 2)
    return res.status(400).json({
      error:
        'The mutation matrix file must include at least one sample column.',
    });
  if (signatureData) {
    const signatureCols = signatureData.length
      ? Object.keys(signatureData[0])
      : [];
    // COSMIC/decomposed signature files use "Type"; de novo files use "MutationType"
    const signatureFirstCol = normalize(signatureCols[0]);
    if (signatureFirstCol !== 'mutationtype' && signatureFirstCol !== 'type')
      return res.status(400).json({
        error: `The signature-profile file's first column must be "MutationType" or "Type" (found "${
          signatureCols[0] ?? ''
        }").`,
      });
    if (signatureCols.length < 2)
      return res.status(400).json({
        error:
          'The signature-profile file must include at least one signature column.',
      });
  }

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
        { signatureSetName, study: 'Reference' },
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
