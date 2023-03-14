import path from 'path';
import logger from '../logger.js';
import formidable from 'formidable';
import fs from 'fs-extra';
import { randomUUID } from 'crypto';
import { validate } from 'uuid';
import Papa from 'papaparse';
import tar from 'tar';
import rWrapper from 'r-wrapper';
import AWS from 'aws-sdk';
import XLSX from 'xlsx';
import replace from 'replace-in-file';
import express from 'express';
import config from '../../config.json' assert { type: 'json' };
import { getObjectBuffer } from '../s3.js';
const r = rWrapper.async;

if (config.aws) AWS.config.update(config.aws);

// config info for R functions
const rConfig = {
  s3Data: config.data.s3,
  bucket: config.data.bucket,
  localData: path.resolve(config.data.localData),
  wd: path.resolve(config.results.folder),
};

function parseCSV(filepath) {
  const file = fs.createReadStream(filepath);
  return new Promise((resolve, reject) => {
    Papa.parse(file, {
      header: true,
      complete(results, file) {
        resolve(results.data);
      },
      error(err, file) {
        logger.info('Error parsing ' + filepath);
        logger.error(err);
        reject(err);
      },
    });
  });
}

async function importUserSession(connection, data, userSchema) {
  const tables = userSchema.filter((e) => !e.type || e.type === 'table');
  const materializedViews = userSchema.filter(
    (e) => e.type === 'materializedView'
  );
  const indexedTables = userSchema.filter((s) => typeof s.index === 'function');
  try {
    // create tables
    for (const { name, schema } of tables) {
      await connection.schema.createTable(name, (table) =>
        schema(table, connection)
      );
    }
    // import data
    for (const [tableName, tableData] of Object.entries(data))
      await connection.batchInsert(tableName, tableData, 100);
    // create "materialized" style tables
    for (const { create } of materializedViews) {
      await create(connection);
    }
    // index tables
    for (const { name, index } of indexedTables) {
      await connection.schema.table(name, index);
    }
    return true;
  } catch (error) {
    logger.error(error);
    return false;
  }
}

function upload(req, res, next) {
  const id = req.params.id;
  if (!validate(id)) next(new Error('Invalid ID'));
  const form = formidable({
    uploadDir: path.resolve(config.folders.input, id),
    multiples: true,
  });

  logger.info(`/upload: Request Project ID:${id}`);
  fs.mkdirSync(form.uploadDir, { recursive: true });

  form
    .on('fileBegin', (field, file) => {
      const destination = path.join(form.uploadDir, file.originalFilename);
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
  const { fn, args, id } = req.body;
  logger.debug('/associationCalc: %o', { ...req.body });
  const sessionId = id || randomUUID();
  // create directory for results if needed
  const savePath = sessionId ? path.join(sessionId, 'results', fn, '/') : null;
  if (sessionId)
    fs.mkdirSync(path.join(rConfig.wd, savePath), { recursive: true });
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

async function getExposureExample(req, res, next) {
  try {
    const { example } = req.params;
    logger.info(`Fetching Exposure example: ${example}`);
    // check exists
    const examplePath = path.resolve(
      config.data.examples,
      'exposure',
      `${example}.tgz`
    );
    if (fs.existsSync(examplePath)) {
      // copy example to results with unique id
      const id = randomUUID();
      const resultsPath = path.resolve(config.results.folder, id);
      await fs.promises.mkdir(resultsPath, { recursive: true });
      // await fs.copy(examplePath, resultsPath);
      await new Promise((resolve, reject) => {
        fs.createReadStream(examplePath)
          .on('end', () => resolve())
          .on('error', (err) => reject(err))
          .pipe(tar.x({ strip: 1, C: resultsPath }));
      });
      const paramsPath = path.join(resultsPath, `params.json`);
      // rename file paths with new ID
      let params = JSON.parse(String(await fs.promises.readFile(paramsPath)));
      const oldID = params.main.id;
      await replace({
        files: paramsPath,
        from: new RegExp(oldID, 'g'),
        to: id,
      });
      params = JSON.parse(String(await fs.promises.readFile(paramsPath)));
      res.json({
        state: {
          ...params,
          main: { ...params.main, id },
        },
      });
    } else {
      throw `Invalid example`;
    }
  } catch (error) {
    next(error);
  }
}

async function getImageS3Batch(req, res, next) {
  // serve static images from s3
  const { keys } = req.body;
  const s3 = new AWS.S3();
  const batch = Object.fromEntries(
    await Promise.all(
      (keys || []).map((Key) => {
        const key = decodeURIComponent(Key);
        return s3
          .getObject({
            Bucket: config.data.bucket,
            Key: key,
          })
          .promise()
          .then(({ Body }) => [
            [path.parse(key).name],
            'data:image/svg+xml;base64,' + Buffer.from(Body).toString('base64'),
          ])
          .catch((error) => {
            logger.error(`${key}: ${error.message}`);
            return [[path.parse(key).name], 'no image available'];
          });
      })
    )
  );
  res.json(batch);
}

async function getImageS3(req, res, next) {
  // serve static images from s3
  const key = req.body.path;
  const s3 = new AWS.S3();
  s3.getObject({
    Bucket: config.data.bucket,
    Key: key,
  })
    .createReadStream()
    .on('error', next)
    .on('pipe', () => res.setHeader('Content-Type', 'image/svg+xml'))
    .pipe(res);
}

async function getFileS3(req, res, next) {
  // serve static files from s3
  const { path } = req.body;
  const s3 = new AWS.S3();
  s3.getObject({
    Bucket: config.data.bucket,
    Key: `msigportal/Database/${path}`,
  })
    .createReadStream()
    .on('error', next)
    .pipe(res);
}

const getDataUsingS3Select = (params) => {
  return new Promise(async (resolve, reject) => {
    try {
      const s3 = new AWS.S3();
      const response = await s3.selectObjectContent(params).promise();
      const events = response.Payload;
      let records = [];
      events.on('data', ({ Records, End }) => {
        if (Records) {
          records.push(Records.Payload);
        } else if (End) {
          let string = Buffer.concat(records).toString('utf8');
          string = string.replace(/\,$/, '');
          const results = JSON.parse(`[${string}]`);
          resolve(results);
        }
      });
      events.on('error', reject);
    } catch (error) {
      reject(error);
    }
  });
};

const router = express.Router();
router.post('/upload/:id?', upload);
router.get('/getExposureExample/:example', getExposureExample);
router.post('/getImageS3Batch', getImageS3Batch);
router.post('/getImageS3', getImageS3);
router.post('/getFileS3', getFileS3);
router.post('/associationWrapper', associationWrapper);

export {
  router,
  parseCSV,
  importUserSession,
  upload,
  getExposureExample,
  getImageS3Batch,
  getImageS3,
  getFileS3,
  associationWrapper,
};
