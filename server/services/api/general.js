import path from 'path';
import formidable from 'formidable';
import fs from 'fs';
import { randomUUID } from 'crypto';
import { validate } from 'uuid';
import Papa from 'papaparse';
import rWrapper from 'r-wrapper';
import Router from 'express-promise-router';
import { mkdirs } from '../utils.js';
const r = rWrapper.async;
import { getObjectBuffer } from '../s3.js';
const env = process.env;

function parseCSV(filepath) {
  const file = fs.createReadStream(filepath);
  return new Promise((resolve, reject) => {
    Papa.parse(file, {
      header: true,
      complete(results, file) {
        resolve(results.data);
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

// async function getImageS3Batch(req, res, next) {
//   // serve static images from s3
//   const { logger } = req.app.locals;
//   const { keys } = req.body;

//   const batch = Object.fromEntries(
//     await Promise.all(
//       keys.map(async (Key) => {
//         const key = decodeURIComponent(Key);
//         try {
//           const image = await getObjectBuffer(key, env.DATA_BUCKET);
//           return [
//             [path.parse(key).name],
//             'data:image/svg+xml;base64,' + image.toString('base64'),
//           ];
//         } catch (error) {
//           logger.error(`${key}: ${error.message}`);
//           logger.error(error);
//           return [
//             [path.parse(key).name],
//             'no image available. ' + error.message,
//           ];
//         }
//       })
//     )
//   );
//   res.json(batch);
// }

// async function getImageS3(req, res, next) {
//   // serve static images from s3
//   const key = req.body?.path;
//   if (!key) {
//     next('Missing path to image');
//   } else {
//     const image = await getObjectBuffer(key, env.DATA_BUCKET);
//     res.setHeader('Content-Type', 'image/svg+xml');
//     res.send(image);
//   }
// }

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

const router = Router();
router.post('/upload/:id?', upload);
// router.post('/getImageS3Batch', getImageS3Batch);
// router.post('/getImageS3', getImageS3);
router.post('/getFileS3', getFileS3);
router.post('/associationWrapper', associationWrapper);

export {
  router,
  parseCSV,
  upload,
  getFileS3,
  associationWrapper,
};
