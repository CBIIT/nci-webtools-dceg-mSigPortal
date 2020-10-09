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
  console.log(svgList);

  return {
    svgList: svgList,
    statistics: statistics,
    matrixList: matrixList,
    downloads: downloads,
  };
}

function getRelativePath(paths) {
  console.log('paths', paths);
  let newPaths = {};
  Object.keys(paths).map((key) => {
    const fullPath = path.resolve(paths[key]);
    newPaths[key] = fullPath.replace(
      path.resolve(config.results.folder) + '/',
      ''
    );
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
}

async function visualizationProfilerExtraction(req, res, next) {
  req.setTimeout(15 * 60 * 1000);
  res.setTimeout(15 * 60 * 1000, () => {
    res.status(504).send('request timed out');
  });

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
    req.body.fn
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
      dataPath: path.join(config.data.folder, 'signature_visualization/'),
    });

    const { stdout, output } = JSON.parse(wrapper);
    // console.log('wrapper return', JSON.parse(wrapper));
    console.log('output', output);
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
      [
        req.body.profileType,
        path.join(config.data.folder, 'signature_visualization/'),
      ]
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
      path.join(config.data.folder, 'signature_visualization/'),
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
      [path.join(config.data.folder, 'signature_visualization/')]
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
    const list = await r('services/R/visualizeWrapper.R', 'getPublicData', [
      req.body.study,
      req.body.cancerType,
      req.body.experimentalStrategy,
      path.join(config.data.folder, 'signature_visualization/'),
    ]);
    logger.info('/getPublicOptions: Complete');

    const projectID = uuidv4();

    res.json({ svgList: JSON.parse(list), projectID: projectID });
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
  form.on('file', (field, file) => {
    const uploadPath = path.join(form.uploadDir, file.name);
    if (field == 'inputFile') form.filePath = uploadPath;
    if (field == 'bedFile') form.bedPath = uploadPath;
    fs.rename(file.path, uploadPath, (err) => {
      if (err) {
        logger.info(`/UPLOAD: Failed to upload file: ${file.name}`);
        logger.error(err);
        res.status(404).json({
          msg: `Failed to upload file: ${file.name}<br><br>
        Review your selected File Type and try again.`,
          error: err,
        });
      } else {
        logger.info(`/UPLOAD: Successfully uploaded file: ${file.name}`);
      }
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
    res.json({
      projectID: projectID,
      filePath: form.filePath,
      bedPath: form.bedPath || '',
    });
  });
}

function getPublicSVG(req, res, next) {
  const svgPath = path.resolve(req.body.path);
  if (svgPath.indexOf(path.resolve(config.data.folder)) == 0) {
    const s = fs.createReadStream(svgPath);

    s.on('open', () => {
      res.set('Content-Type', 'image/svg+xml');
      s.pipe(res);
      logger.debug(`/getPublicSVG: Serving ${svgPath}`);
    });
    s.on('error', () => {
      res.set('Content-Type', 'text/plain');
      res.status(500).end('Not found');
      logger.debug(`/getPublicSVG: Serving Error retrieving ${svgPath}`);
    });
  } else {
    logger.info('traversal error');
    res.status(500).end('Not found');
  }
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
  const projectID = req.body.projectID || uuidv4();
  const savePath = path.join(
    config.results.folder,
    projectID,
    'results',
    req.body.fn
  );

  fs.mkdirSync(savePath, { recursive: true });

  const wrapper = await r('services/R/exploringWrapper.R', req.body.fn, {
    ...req.body.args,
    projectID: projectID,
    pythonOutput: path.join(
      config.results.folder,
      req.body.projectID,
      'results/output'
    ),
    savePath: savePath,
    dataPath: path.join(config.data.folder, 'signature_visualization/'),
  }).catch(next);

  const { stdout, output } = JSON.parse(wrapper);

  return res.json({
    debugR: stdout,
    output: getRelativePath(output),
    projectID: projectID,
  });
}

async function getReferenceSignatureData(req, res, next) {
  logger.info('/getReferenceSignatures: Request');

  const data = await r(
    'services/R/exploringWrapper.R',
    'getReferenceSignatureData',
    {
      ...req.body,
      dataPath: path.join(config.data.folder, 'signature_visualization/'),
    }
  ).catch(next);

  logger.info('/getReferenceSignatures: Success');
  return res.json(JSON.parse(data));
}

async function submitQueue(req, res, next) {
  const projectID = req.body.args.projectID[1];
  const date = new Date();
  const isoDate = new Date(
    date.getTime() - date.getTimezoneOffset() * 60000
  ).toISOString();
  const day = isoDate.split('T')[0];
  const time = isoDate.split('T')[1].split('Z')[0].substring(0, 5);
  const sqs = new AWS.SQS();

  try {
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
          timestamp: `${day} ${time} UTC`,
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

      logger.info('/fetchResults: Found Params')
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
  getPublicSVG,
  download,
  exploringR,
  getReferenceSignatureData,
  submitQueue,
  fetchResults,
};
