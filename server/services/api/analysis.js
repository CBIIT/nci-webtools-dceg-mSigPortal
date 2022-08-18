const path = require('path');
const logger = require('../logger');
const formidable = require('formidable');
const fs = require('fs-extra');
const { randomUUID } = require('crypto');
const Papa = require('papaparse');
const tar = require('tar');
const r = require('r-wrapper').async;
const AWS = require('aws-sdk');
const XLSX = require('xlsx');
const replace = require('replace-in-file');
const express = require('express');
const config = require('../../config.json');

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
  const tables = userSchema.filter(({ type }) => !type || type === 'table');
  const materializedViews = userSchema.filter(
    ({ type }) => type === 'materializedView'
  );
  const indexedTables = userSchema.filter((s) => typeof s.index === 'function');
  const tableName = `seqmatrix`;

  try {
    for (const { name, schema } of tables) {
      await connection.schema.createTable(name, (table) =>
        schema(table, connection)
      );
    }

    await connection.batchInsert(tableName, data, 100);

    for (const { create } of materializedViews) {
      await create(connection);
    }

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
  const projectID = randomUUID();
  const form = formidable({
    uploadDir: path.join(config.results.folder, projectID),
    multiples: true,
  });

  logger.info(`/upload: Request Project ID:${projectID}`);

  fs.mkdirSync(form.uploadDir);

  form
    .on('fileBegin', (field, file) => {
      uploadPath = path.join(form.uploadDir, file.name);
      if (field == 'inputFile') form.filePath = uploadPath;
      if (field == 'bedFile') form.bedPath = uploadPath;
      if (field == 'exposureFile') form.exposurePath = uploadPath;

      file.path = uploadPath;
    })
    .on('error', (err) => {
      logger.info('/UPLOAD: An error occured\n' + err);
      logger.error(err);
      res.status(500).json({
        msg: 'An error occured while trying to upload',
        err: err,
      });
    })
    .on('end', () => {
      res.json({
        projectID: projectID,
        filePath: form.filePath,
        bedPath: form.bedPath,
        exposurePath: form.exposurePath,
      });
    });

  form.parse(req);
}

async function explorationWrapper(req, res, next) {
  const { fn, args, projectID: id, type = 'calc' } = req.body;
  logger.debug('/explorationWrapper: %o', { ...req.body });

  const projectID = id ? id : type == 'calc' ? randomUUID() : false;
  // create directory for results if needed
  const savePath = projectID ? path.join(projectID, 'results', fn, '/') : null;
  if (projectID)
    fs.mkdirSync(path.join(rConfig.wd, savePath), { recursive: true });

  try {
    const wrapper = await r('services/R/explorationWrapper.R', 'wrapper', {
      fn,
      args,
      config: {
        ...rConfig,
        savePath,
        projectID,
      },
    });

    const { stdout, ...rest } = JSON.parse(wrapper);

    logger.debug(stdout);

    res.json({
      projectID,
      stdout,
      ...rest,
    });
  } catch (err) {
    logger.error(`/explorationCalc: An error occured with fn: ${fn}`);
    next(err);
  }
}

async function associationWrapper(req, res, next) {
  const { fn, args, projectID: id } = req.body;
  logger.debug('/associationCalc: %o', { ...req.body });

  const projectID = id ? id : randomUUID();

  // create directory for results if needed
  const savePath = projectID ? path.join(projectID, 'results', fn, '/') : null;
  if (projectID)
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
      projectID,
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
      const oldID = params.main.projectID;
      await replace({
        files: paramsPath,
        from: new RegExp(oldID, 'g'),
        to: id,
      });
      params = JSON.parse(String(await fs.promises.readFile(paramsPath)));

      res.json({
        state: {
          ...params,
          main: { ...params.main, projectID: id },
        },
      });
    } else {
      throw `Invalid example`;
    }
  } catch (error) {
    next(error);
  }
}

// Publications page data
async function getPublications(req, res, next) {
  let buffers = [];
  const filestream = new AWS.S3()
    .getObject({
      Bucket: config.data.bucket,
      Key: `${config.data.s3}Others/Publications.xlsx`,
    })
    .createReadStream();

  filestream
    .on('data', (data) => buffers.push(data))
    .on('end', () => {
      const buffer = Buffer.concat(buffers);
      const workbook = XLSX.read(buffer);
      excelToJSON(workbook);
    })
    .on('error', next);

  function excelToJSON(workbook) {
    const sheetNames = workbook.SheetNames;
    const data = sheetNames.reduce(
      (acc, sheet) => ({
        ...acc,
        [sheet]: XLSX.utils.sheet_to_json(workbook.Sheets[sheet]),
      }),
      {}
    );

    res.json(data);
  }
}

async function getImageS3Batch(req, res, next) {
  // serve static images from s3
  const { keys } = req.body;
  const s3 = new AWS.S3();

  const batch = await Promise.all(
    (keys || []).map((Key) =>
      s3
        .getObject({
          Bucket: config.data.bucket,
          Key,
        })
        .promise()
        .then(
          ({ Body }) =>
            'data:image/svg+xml;base64,' + Buffer.from(Body).toString('base64')
        )
        .catch((err) => '')
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

router.post('/upload', upload);

router.post('/explorationWrapper', explorationWrapper);

router.get('/getExposureExample/:example', getExposureExample);

router.get('/getPublications', getPublications);

router.post('/getImageS3Batch', getImageS3Batch);

router.post('/getImageS3', getImageS3);

router.post('/getFileS3', getFileS3);

router.post('/associationWrapper', associationWrapper);

module.exports = {
  router,
  parseCSV,
  importUserSession,
  upload,
  explorationWrapper,
  getExposureExample,
  getPublications,
  getImageS3Batch,
  getImageS3,
  getFileS3,
  associationWrapper,
};
