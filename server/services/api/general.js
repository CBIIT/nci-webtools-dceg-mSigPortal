import path from 'path';
import formidable from 'formidable';
import fs from 'fs';
import { randomUUID } from 'crypto';
import { validate } from 'uuid';
import Papa from 'papaparse';
import rWrapper from 'r-wrapper';
import Router from 'express-promise-router';
import archiver from 'archiver';
import { mkdirs } from '../utils.js';
const r = rWrapper.async;
import { getObjectBuffer } from '../s3.js';
const env = process.env;

function parseCSV(filepath) {
  const file = fs.createReadStream(filepath);
  return new Promise((resolve, reject) => {
    Papa.parse(file, {
      header: true,
      skipEmptyLines: true,
      complete(results, file) {
        resolve(results.data);
      },
      error(err, file) {
        reject(err);
      },
    });
  });
}

export function parseTSV(filepath) {
  const file = fs.createReadStream(filepath);
  return new Promise((resolve, reject) => {
    Papa.parse(file, {
      header: true,
      skipEmptyLines: true,
      delimiter: '\t',
      complete(results, file) {
        resolve(results);
      },
      error(err, file) {
        reject(err);
      },
    });
  });
}

function upload(req, res, next) {
  const { logger } = req.app.locals;
  const { id } = req.params;
  if (!validate(id)) next(new Error('Invalid ID'));

  const form = formidable({
    uploadDir: path.resolve(env.INPUT_FOLDER, id),
    multiples: true,
  });

  logger.info(`/upload: Request Project ID:${id}`);
  fs.mkdirSync(form.uploadDir, { recursive: true });

  form
    .on('fileBegin', (field, file) => {
      const destination = path.resolve(form.uploadDir, file.originalFilename);
      file.filepath = destination;
    })
    .on('error', (err) => {
      logger.info('/UPLOAD: An error occurred\n' + err);
      logger.error(err);
      res.status(500).json({
        msg: 'An error occurred during upload',
        err: err,
      });
    })
    .on('end', () => {
      res.json({ id });
    });

  form.parse(req);
}

async function associationWrapper(req, res, next) {
  const { logger } = req.app.locals;
  const { fn, args, id } = req.body;
  const sessionId = id || randomUUID();

  // config info for R functions
  const rConfig = {
    prefix: env.DATA_BUCKET_PREFIX,
    bucket: env.DATA_BUCKET,
    wd: path.resolve(env.DATA_FOLDER),
  };

  // create directory for results if needed
  const savePath = sessionId
    ? path.join('output', sessionId, 'results', fn, '/')
    : '';
  if (sessionId) await mkdirs([path.join(rConfig.wd, savePath)]);

  try {
    const wrapper = await r('services/R/associationWrapper.R', 'wrapper', {
      fn,
      args,
      config: {
        ...rConfig,
        savePath: savePath,
      },
    });
    const { stdout, ...rest } = JSON.parse(wrapper);
    logger.debug(stdout);
    res.json({
      sessionId,
      stdout,
      ...rest,
    });
  } catch (err) {
    logger.error(`/associationWrapper: An error occured with fn:${fn}`);
    next(err);
  }
}

async function getFileS3(req, res, next) {
  // serve static files from s3
  const { path } = req.body;
  if (path) {
    const file = await getObjectBuffer(
      `msigportal/Database/${path}`,
      env.DATA_BUCKET
    );
    res.send(file);
  } else {
    next('Missing path to file');
  }
}

export async function downloadOutput(req, res, next) {
  const { id } = req.params;
  const output = path.resolve(env.OUTPUT_FOLDER, id);
  const archive = archiver('zip', { zlib: { level: 6 } });

  if (!validate(id)) res.status(500).json(`${id} is not a valid ID`);
  if (!fs.existsSync(output)) res.status(500).json(`${id} does not exist`);

  res.attachment(`${id}.zip`);
  archive.directory(output, false).pipe(res);
  archive.finalize();
}

const router = Router();
router.post('/upload/:id?', upload);
router.get('/downloadOutput/:id', downloadOutput);
router.post('/getFileS3', getFileS3);
router.post('/associationWrapper', associationWrapper);

export { router, parseCSV };
