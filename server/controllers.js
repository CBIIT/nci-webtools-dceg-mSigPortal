const path = require('path');
const logger = require('./logger');
const { spawn } = require('promisify-child-process');
const formidable = require('formidable');
const fs = require('fs');
const { v4: uuidv4, validate } = require('uuid');
const Papa = require('papaparse');
const tar = require('tar');
const r = require('r-wrapper').async;
const AWS = require('aws-sdk');
const config = require('./config.json');

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

async function getSummaryFiles(resultsPath) {
  const svgListPath = path.join(resultsPath, 'svg_files_list.txt');
  const statisticsPath = path.join(resultsPath, 'Statistics.txt');
  const matrixPath = path.join(resultsPath, 'matrix_files_list.txt');
  const downloadsPath = path.join(resultsPath, 'output');
  let matrixList = [];
  let statistics = '';
  let downloads = [];
  let svgList = await parseCSV(svgListPath);

  if (fs.existsSync(matrixPath)) matrixList = await parseCSV(matrixPath);
  if (fs.existsSync(statisticsPath))
    statistics = fs.readFileSync(statisticsPath, 'utf8');

  if (fs.existsSync(downloadsPath)) {
    downloads = fs
      .readdirSync(downloadsPath)
      .filter((file) => file.endsWith('.tar.gz'))
      .map((file) => file);
  }

  svgList.forEach(
    (plot) => (plot.Path = getRelativePath({ Path: plot.Path }).Path)
  );

  return {
    svgList: svgList,
    statistics: statistics,
    matrixList: matrixList,
    downloads: downloads,
  };
}

function getRelativePath(paths) {
  let newPaths = {};
  const resultsPath = path.resolve(config.results.folder);
  const dataPath = path.resolve(config.data.database);

  Object.keys(paths).map((key) => {
    const fullPath = path.resolve(paths[key]);
    if (fullPath.includes(resultsPath))
      newPaths[key] = fullPath.replace(resultsPath + '/', '');
    if (fullPath.includes(dataPath))
      newPaths[key] = fullPath.replace(dataPath + '/', '');
  });
  return newPaths;
}

async function profilerExtraction(params) {
  logger.info('/profilerExtraction: Spawning Python Process');
  // update path
  params.outputDir[1] = path.join(
    config.results.folder,
    params.outputDir[1],
    'results'
  );

  const args = Object.values(params);
  const cli = args.reduce((params, arg) => [...params, ...arg]);

  logger.debug('/profilerExtraction: CLI args\n' + cli);

  try {
    const { stdout, stderr } = await spawn(
      'python3',
      ['services/python/mSigPortal_Profiler_Extraction.py', ...cli],
      { encoding: 'utf8' }
    );

    return {
      stdout: stdout,
      stderr: stderr,
      projectPath: path.join(config.results.folder, params.projectID[1]),
    };
  } catch (e) {
    throw e;
  }
}

async function visualizationProfilerExtraction(req, res, next) {
  req.setTimeout(15 * 60 * 1000);
  res.setTimeout(15 * 60 * 1000, () => {
    res.status(504).send('request timed out');
  });

  try {
    const { stdout, stderr, projectPath } = await profilerExtraction(req.body);
    const resultsPath = path.join(projectPath, 'results');

    if (fs.existsSync(path.join(resultsPath, 'svg_files_list.txt'))) {
      res.json({ stdout, stderr, ...(await getSummaryFiles(resultsPath)) });
    } else {
      logger.info(
        '/profilerExtraction: An Error Occured While Extracting Profiles'
      );
      res.status(500).json({
        msg:
          'An error occured durring profile extraction. Please review your input parameters and try again.',
        stdout,
        stderr,
      });
    }
  } catch (err) {
    next(err);
  }
}

async function getSummary(req, res, next) {
  logger.info('/getSummary: Retrieving Summary');
  console.log('summary', req.body);
  const resultsPath = path.join(
    config.results.folder,
    req.body.projectID,
    'results'
  );

  if (fs.existsSync(path.join(resultsPath, 'svg_files_list.txt'))) {
    res.json(await getSummaryFiles(resultsPath));
  } else {
    logger.info('/getSummary: Summary file not found');
    res.status(500).json({
      msg: 'Summary file not found',
    });
  }
}

async function visualizeR(req, res, next) {
  logger.info('/visualizeR: function ' + req.body.fn);
  console.log('args', req.body);
  const savePath = path.join(
    config.results.folder,
    req.body.projectID,
    'results',
    req.body.fn,
    '/'
  );

  fs.mkdirSync(savePath, { recursive: true });

  try {
    const wrapper = await r('services/R/visualizeWrapper.R', req.body.fn, {
      ...req.body.args,
      projectID: req.body.projectID,
      pythonOutput: path.join(
        config.results.folder,
        req.body.projectID,
        'results/output'
      ),
      savePath: savePath,
      dataPath: path.join(config.data.database),
    });

    const { stdout, output } = JSON.parse(wrapper);
    // logger.debug(stdout);

    res.json({
      debugR: stdout,
      output: getRelativePath(output),
    });
  } catch (err) {
    logger.info('/visualizeR: An error occured');
    logger.error(err);
    res.status(500).json(err.message);
  }
}

async function getReferenceSignatureSets(req, res, next) {
  logger.info('/getReferenceSignatureSets: Calling R Wrapper');
  console.log('args', req.body);

  try {
    const list = await r(
      'services/R/visualizeWrapper.R',
      'getReferenceSignatureSets',
      [req.body.profileType, path.join(config.data.database)]
    );

    // console.log('SignatureReferenceSets', list);

    res.json(list);
  } catch (err) {
    logger.info('/getReferenceSignatureSets: An error occured');
    logger.error(err);
    res.status(500).json(err.message);
  }
}

async function getSignatures(req, res, next) {
  logger.info('/getSignatures: Calling R Wrapper');

  try {
    const list = await r('services/R/visualizeWrapper.R', 'getSignatures', [
      req.body.profileType,
      req.body.signatureSetName,
      path.join(config.data.database),
    ]);

    // console.log('signatures', list);

    res.json(list);
  } catch (err) {
    logger.info('/getSignatures: An error occured');
    logger.error(err);
    res.status(500).json(err.message);
  }
}

async function getPublicDataOptions(req, res, next) {
  logger.info('/getPublicOptions: Request');
  try {
    const list = await r(
      'services/R/visualizeWrapper.R',
      'getPublicDataOptions',
      [path.join(config.data.database)]
    );

    res.json(JSON.parse(list));
    logger.info('/getPublicOptions: Success');
  } catch (err) {
    logger.info('/getPublicOptions: An error occured');
    logger.error(err);
    res.status(500).json(err.message);
  }
}

async function getPublicData(req, res, next) {
  logger.info('/getPublicOptions: Calling R Wrapper');
  try {
    const projectID = uuidv4();
    const list = await r('services/R/visualizeWrapper.R', 'getPublicData', [
      req.body.study,
      req.body.cancerType,
      req.body.experimentalStrategy,
      path.join(config.data.database),
    ]);
    logger.info('/getPublicOptions: Complete');

    let svgList = JSON.parse(list);
    svgList.forEach(
      (plot) => (plot.Path = getRelativePath({ Path: plot.Path }).Path)
    );

    res.json({
      svgList: svgList,
      projectID: projectID,
    });
  } catch (err) {
    logger.info('/getPublicOptions: An error occured');
    logger.error(err);
    res.status(500).json(err.message);
  }
}

function upload(req, res, next) {
  const projectID = uuidv4();
  const form = formidable({
    uploadDir: path.join(config.results.folder, projectID),
    multiples: true,
  });

  logger.info(`/upload: Request Project ID:${projectID}`);

  fs.mkdirSync(form.uploadDir);

  form.parse(req);
  form.on('file', async (field, file) => {
    const uploadPath = path.join(form.uploadDir, file.name);
    if (field == 'inputFile') form.filePath = uploadPath;
    if (field == 'bedFile') form.bedPath = uploadPath;

    let data = await fs.promises.readFile(file.path);
    await fs.promises
      .writeFile(uploadPath, data)
      .then(logger.info(`/UPLOAD: Successfully uploaded file: ${file.name}`))
      .catch((err) => {
        logger.info(`/UPLOAD: Failed to upload file: ${file.name}`);
        logger.error(err);
        res.status(404).json({
          msg: `Failed to upload file: ${file.name}<br><br>
        Review your selected File Type and try again.`,
          error: err,
        });
      });
  });
  form.on('error', (err) => {
    logger.info('/UPLOAD: An error occured\n' + err);
    logger.error(err);
    res.status(500).json({
      msg: 'An error occured while trying to upload',
      err: err,
    });
  });
  form.on('end', () => {
    logger.info('Upload Complete');
    res.json({
      projectID: projectID,
      filePath: form.filePath,
      bedPath: form.bedPath || '',
    });
  });
}

function download(req, res, next) {
  logger.info(`/visualize/download: id:${req.query.id} file:${req.query.file}`);
  const file = path.resolve(
    path.join(
      config.results.folder,
      req.query.id,
      'results/output',
      req.query.file
    )
  );
  if (file.indexOf(path.resolve(config.results.folder)) == 0) {
    res.download(file);
  } else {
    logger.info('traversal error');
    res.status(500).end('Not found');
  }
}

async function exploringR(req, res, next) {
  logger.info('/exploringR: function ' + req.body.fn);
  console.log('args', req.body);
  const projectID = req.body.projectID ? req.body.projectID : uuidv4();
  const rootDir = path.join(config.results.folder, projectID);
  const savePath = path.join(rootDir, 'results', req.body.fn, '/');

  fs.mkdirSync(savePath, { recursive: true });

  try {
    const wrapper = await r('services/R/exploringWrapper.R', req.body.fn, {
      ...req.body.args,
      projectID: projectID,
      pythonOutput: path.join(
        config.results.folder,
        req.body.projectID,
        'results/output'
      ),
      rootDir: rootDir,
      savePath: savePath,
      dataPath: path.join(config.data.database),
    });

    const { stdout, output } = JSON.parse(wrapper);

    res.json({
      debugR: stdout,
      output: getRelativePath(output),
      projectID: projectID,
    });
  } catch (err) {
    logger.info(req.body.fn + ' failed');
    res.json({ debugR: err.stderr });
  }
}

async function getReferenceSignatureData(req, res, next) {
  logger.info('/getReferenceSignatures: Request');

  const data = await r(
    'services/R/exploringWrapper.R',
    'getReferenceSignatureData',
    {
      ...req.body,
      dataPath: path.join(config.data.database),
    }
  ).catch(next);

  logger.info('/getReferenceSignatures: Success');
  res.json(JSON.parse(data));
}

async function submitQueue(req, res, next) {
  const projectID = req.body.args.projectID[1];
  const sqs = new AWS.SQS();

  try {
    // upload archived project directory
    await new AWS.S3()
      .upload({
        Body: tar.c({ gzip: true, C: config.results.folder }, [projectID]),
        Bucket: config.s3.bucket,
        Key: `${config.s3.outputKeyPrefix}${projectID}/${projectID}.tgz`,
      })
      .promise();

    const { QueueUrl } = await sqs
      .getQueueUrl({ QueueName: config.queue.url })
      .promise();

    await sqs
      .sendMessage({
        QueueUrl: QueueUrl,
        MessageDeduplicationId: projectID,
        MessageGroupId: projectID,
        MessageBody: JSON.stringify({
          ...req.body,
          timestamp: new Date().toLocaleString('en-US', {
            timeZone: 'America/New_York',
          }),
        }),
      })
      .promise();

    logger.info('Queue submitted ID: ' + projectID);
    res.json({ projectID });
  } catch (err) {
    logger.info('Queue failed to submit ID: ' + projectID);
    next(err);
  }
}

async function fetchResults(req, res, next) {
  try {
    const s3 = new AWS.S3();
    const { id } = req.params;

    // validate id format
    if (!validate(id)) {
      throw `Invalid id`;
    }

    // ensure output directory exists
    const resultsFolder = path.resolve(config.results.folder, id);
    await fs.promises.mkdir(resultsFolder, { recursive: true });

    // find objects which use the specified id as the prefix
    const objects = await s3
      .listObjectsV2({
        Bucket: config.s3.bucket,
        Prefix: `${config.s3.outputKeyPrefix}${id}/`,
      })
      .promise();

    // download results
    for (let { Key } of objects.Contents) {
      const filename = path.basename(Key);
      const filepath = path.resolve(resultsFolder, filename);

      // download results if they do not exist
      if (!fs.existsSync(filepath)) {
        logger.info(`Downloading result: ${Key}`);
        const object = await s3
          .getObject({
            Bucket: config.s3.bucket,
            Key,
          })
          .promise();

        await fs.promises.writeFile(filepath, object.Body);
        // extract and delete archive
        if (path.extname(filename) == '.tgz') {
          fs.createReadStream(filepath)
            .pipe(tar.x({ strip: 1, C: resultsFolder }))
            .once('finish', () => fs.unlink(filepath, next));
        }
      }
    }

    let paramsFilePath = path.resolve(resultsFolder, `params.json`);

    if (fs.existsSync(paramsFilePath)) {
      const params = JSON.parse(
        String(await fs.promises.readFile(paramsFilePath))
      );

      logger.info('/fetchResults: Found Params');
      res.json(params);
    } else {
      throw `Invalid id`;
    }
  } catch (error) {
    next(error);
  }
}

module.exports = {
  profilerExtraction,
  visualizationProfilerExtraction,
  getSummary,
  visualizeR,
  getReferenceSignatureSets,
  getSignatures,
  getPublicDataOptions,
  getPublicData,
  upload,
  download,
  exploringR,
  getReferenceSignatureData,
  submitQueue,
  fetchResults,
  parseCSV,
};
