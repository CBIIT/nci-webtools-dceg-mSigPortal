const path = require('path');
const logger = require('./logger');
const { spawn } = require('promisify-child-process');
const formidable = require('formidable');
const fs = require('fs-extra');
const { v4: uuidv4, validate } = require('uuid');
const Papa = require('papaparse');
const tar = require('tar');
const r = require('r-wrapper').async;
const AWS = require('aws-sdk');
const XLSX = require('xlsx');
const replace = require('replace-in-file');
const config = require('./config.json');

if (config.aws) AWS.config.update(config.aws);

const dataArgs = {
  s3Data: config.data.s3,
  localData: path.join(config.data.localData),
  bucket: config.data.bucket,
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

// convert array of objects into 2d array
function to2dArray(list) {
  if (list)
    return {
      columns: Object.keys(list[0]),
      data: list.map((obj) => Object.values(obj)),
    };
}

// retrieves parsed data files - modify paths if needed
async function getResultsFiles(resultsPath, id = '') {
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
    (plot) => (plot.Path = getRelativePath({ Path: plot.Path }, id).Path)
  );

  return {
    svgList: to2dArray(svgList),
    statistics: statistics,
    matrixList: to2dArray(matrixList),
    downloads: downloads,
  };
}

function getRelativePath(paths, id = '') {
  let newPaths = {};
  const resultsPath = path.resolve(config.results.folder, id);
  const dataPath = path.resolve(config.data.localData);

  Object.keys(paths).map((key) => {
    const fullPath = path.resolve(paths[key]);

    if (fullPath.includes(resultsPath))
      newPaths[key] = fullPath.replace(resultsPath, '');
    else if (fullPath.includes(dataPath))
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
  } catch (error) {
    // const error = { code, signal, stdout, stderr };
    logger.info('Profiler Extraction Error');
    throw error;
  }
}

async function visualizationProfilerExtraction(req, res, next) {
  res.setTimeout(15 * 60 * 1000, () => {
    res.status(504).send('request timed out');
  });

  try {
    const { stdout, stderr, projectPath } = await profilerExtraction(req.body);
    const resultsPath = path.join(projectPath, 'results');

    if (fs.existsSync(path.join(resultsPath, 'svg_files_list.txt'))) {
      res.json({
        stdout,
        stderr,
        ...(await getResultsFiles(resultsPath, req.body.projectID[1])),
      });
    } else {
      logger.info(
        '/profilerExtraction: An Error Occured While Extracting Profiles'
      );
      res.status(500).json({
        stdout,
        stderr,
      });
    }
  } catch ({ stdout, stderr }) {
    res.status(500).json({
      stdout,
      stderr,
    });
  }
}

async function getResults(req, res, next) {
  logger.info(`/getResults: Retrieving Results for ${req.body.projectID}`);

  const userResults = path.resolve(
    config.results.folder,
    req.body.projectID,
    'results'
  );

  if (fs.existsSync(path.join(userResults, 'svg_files_list.txt'))) {
    res.json(await getResultsFiles(userResults, req.body.projectID));
  } else {
    logger.info('/getResults: Results not found');
    res.status(500).json('Results not found');
  }
}

// Visualization Calculation functions
async function visualizationCalc(req, res, next) {
  const { fn, args, projectID } = req.body;
  logger.debug(`/visualizationCalc: %o`, req.body);

  // create save directory if needed
  const savePath = path.join(
    config.results.folder,
    projectID,
    'results',
    fn,
    '/'
  );
  fs.mkdirSync(savePath, { recursive: true });

  try {
    const wrapper = await r('services/R/visualizeWrapper.R', fn, {
      ...args,
      projectID: projectID,
      pythonOutput: path.join(
        config.results.folder,
        projectID,
        'results/output'
      ),
      savePath: savePath,
      ...dataArgs,
    });

    const { stdout, output, ...rest } = JSON.parse(wrapper);
    // logger.debug(stdout);

    res.json({
      ...rest,
      debugR: stdout,
      output: getRelativePath(output, projectID),
    });
  } catch (err) {
    logger.info(`/visualizationCalc: An error occured with fn:${fn}`);
    res.status(500).json(err.message);
    next(err);
  }
}

// Visualization data/util functions
async function visualizationData(req, res, next) {
  const { fn, args } = req.body;
  logger.debug(`/visualizationData: %o`, req.body);

  try {
    const data = await r('services/R/visualizeWrapper.R', fn, {
      ...args,
      ...dataArgs,
    });

    res.json(data);
  } catch (err) {
    logger.info(`/visualizationData: An error occured with fn:${fn}`);
    res.status(500).json(err.message);
    next(err);
  }
}

async function getSignaturesUser(req, res, next) {
  logger.info('/getSignaturesUser: Parsing File');
  try {
    const file = path.resolve(req.body.path);
    logger.debug(file);
    if (file.indexOf(path.resolve(config.results.folder)) == 0) {
      const data = await parseCSV(file);
      res.json(to2dArray(data));
    } else {
      logger.info('traversal error');
      res.status(500).end('Not found');
    }
  } catch (err) {
    logger.info('/getSignaturesUser: An error occured');
    logger.error(err);
    res.status(500).json(err.message);
  }
}

async function getPublicData(req, res, next) {
  logger.info('/getPublicOptions: Calling R Wrapper');
  try {
    const projectID = uuidv4();
    const list = await r('services/R/visualizeWrapper.R', 'getPublicData', {
      study: req.body.study,
      cancerType: req.body.cancerType,
      experimentalStrategy: req.body.experimentalStrategy,
      ...dataArgs,
    });
    logger.info('/getPublicOptions: Complete');

    let svgList = JSON.parse(list);

    res.json({
      svgList: to2dArray(svgList),
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

function visualizationDownload(req, res, next) {
  logger.info(
    `/visualization/download: id:${req.query.id} file:${req.query.file}`
  );
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

// Generate public data files for download
async function visualizationDownloadPublic(req, res, next) {
  logger.info(`/visualization/downloadPublic`);
  const { id, ...args } = req.body.args;
  const { study, cancerType, experimentalStrategy } = args;
  const savePath = path.join(
    config.results.folder,
    id,
    `results/msigportal-${study}-${cancerType}-${experimentalStrategy}`
  );

  try {
    await r('services/R/visualizeWrapper.R', 'downloadPublicData', {
      ...args,
      ...dataArgs,
      savePath,
    });

    const file = path.resolve(
      path.join(
        config.results.folder,
        id,
        `results/msigportal-${study}-${cancerType}-${experimentalStrategy}.tar.gz`
      )
    );

    if (file.indexOf(path.resolve(config.results.folder)) == 0) {
      res.download(file);
    } else {
      logger.error('visualizationDownloadPublic failed');
      res.status(500).end('Not found');
    }
  } catch (err) {
    logger.error(err);
    res.status(500).json(err.message);
    next(err);
  }
}

async function explorationCalc(req, res, next) {
  const { fn, args, projectID = uuidv4() } = req.body;
  logger.debug('/explorationCalc: %o', { ...req.body, projectID });

  const rootDir = path.join(config.results.folder, projectID);
  // create save directory if needed
  const savePath = path.join(rootDir, 'results', fn, '/');
  fs.mkdirSync(savePath, { recursive: true });

  try {
    const wrapper = await r('services/R/explorationWrapper.R', fn, {
      ...args,
      projectID: projectID,
      pythonOutput: path.join(
        config.results.folder,
        projectID,
        'results/output'
      ),
      rootDir: rootDir,
      savePath: savePath,
      ...dataArgs,
    });

    const { stdout, output, ...rest } = JSON.parse(wrapper);

    res.json({
      ...rest,
      debugR: stdout,
      output: { ...output, ...getRelativePath(output, projectID) },
      projectID: projectID,
    });
  } catch (err) {
    logger.info(`/explorationCalc: An error occured with fn:${fn}`);
    res.json({ debugR: err.stderr });
    next(err);
  }
}

async function explorationData(req, res, next) {
  const { fn, args } = req.body;
  logger.debug('/explorationData: %o', req.body);

  try {
    const data = await r('services/R/explorationWrapper.R', fn, {
      args: args,
      ...dataArgs,
    });

    res.json(JSON.parse(data));
  } catch (err) {
    logger.info(`/explorationData: An error occured with fn:${fn}`);
    res.json({ debugR: err.stderr });
    next(err);
  }
}

async function submitQueue(req, res, next) {
  const projectID = req.body.args.projectID[1];
  const sqs = new AWS.SQS();

  try {
    // upload archived project directory
    await new AWS.S3()
      .upload({
        Body: tar
          .c({ sync: true, gzip: true, C: config.results.folder }, [projectID])
          .read(),
        Bucket: config.queue.bucket,
        Key: `${config.queue.inputKeyPrefix}${projectID}/${projectID}.tgz`,
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

async function getQueueResults(req, res, next) {
  try {
    const s3 = new AWS.S3();
    const { id } = req.params;

    logger.info(`Fetch Queue Result: ${id}`);

    // validate id format
    if (!validate(id)) next(new Error(`Invalid request`));

    // ensure output directory exists
    const resultsFolder = path.resolve(config.results.folder, id);
    await fs.promises.mkdir(resultsFolder, { recursive: true });

    // find objects which use the specified id as the prefix
    const objects = await s3
      .listObjectsV2({
        Bucket: config.queue.bucket,
        Prefix: `${config.queue.outputKeyPrefix}${id}/`,
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
            Bucket: config.queue.bucket,
            Key,
          })
          .promise();

        await fs.promises.writeFile(filepath, object.Body);
        // extract and delete archive
        if (path.extname(filename) == '.tgz') {
          await new Promise((resolve, reject) => {
            fs.createReadStream(filepath)
              .on('end', () =>
                fs.unlink(filepath, (err) => {
                  if (err) {
                    reject(err);
                  } else {
                    resolve();
                  }
                })
              )
              .pipe(tar.x({ strip: 1, C: resultsFolder }));
          });
        }
      }
    }

    let paramsPath = path.resolve(resultsFolder, `params.json`);

    if (fs.existsSync(paramsPath)) {
      const data = JSON.parse(String(await fs.promises.readFile(paramsPath)));

      res.json(data);
    } else {
      next(new Error(`Params not found`));
    }
  } catch (error) {
    next(error);
  }
}

async function getVisExample(req, res, next) {
  try {
    const { example } = req.params;
    logger.info(`Fetching example: ${example}`);

    // check exists
    const examplePath = path.resolve(config.data.examples, example);
    const paramsPath = path.join(examplePath, `params.json`);

    if (fs.existsSync(paramsPath)) {
      const params = JSON.parse(String(await fs.promises.readFile(paramsPath)));

      // copy example to results with unique id
      const id = uuidv4();
      const resultsPath = path.resolve(config.results.folder, id);
      await fs.promises.mkdir(resultsPath, { recursive: true });
      await fs.copy(examplePath, resultsPath);

      // rename file paths with new ID if needed
      const svgPath = path.join(resultsPath, 'results', 'svg_files_list.txt');
      if (fs.existsSync(svgPath)) {
        const oldID = params.args.projectID[1];
        const matrixPath = path.join(
          resultsPath,
          'results',
          'matrix_files_list.txt'
        );

        await replace({
          files: [svgPath, matrixPath],
          from: new RegExp(`/${oldID}/`, 'g'),
          to: `/${id}/`,
        });
      }

      res.json({ projectID: id, state: params.state });
    } else {
      throw `Invalid example`;
    }
  } catch (error) {
    next(error);
  }
}

async function getExposureExample(req, res, next) {
  try {
    const { example } = req.params;
    logger.info(`Fetching Exposure example: ${example}`);

    // check exists
    const examplePath = path.resolve(config.data.examples, example);
    const paramsPath = path.join(examplePath, `params.json`);
    console.log(paramsPath);

    if (fs.existsSync(paramsPath)) {
      const params = JSON.parse(String(await fs.promises.readFile(paramsPath)));

      // copy example to results with unique id
      const id = uuidv4();
      const resultsPath = path.resolve(config.results.folder, id);
      await fs.promises.mkdir(resultsPath, { recursive: true });
      await fs.copy(examplePath, resultsPath);

      res.json({ projectID: id, state: params.state });
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

async function getImageS3(req, res, next) {
  // serve static images from s3
  const key = req.body.path;
  const s3 = new AWS.S3();

  res.setHeader('Content-Type', 'image/svg+xml');
  s3.getObject({
    Bucket: config.data.bucket,
    Key: key,
  })
    .createReadStream()
    .on('error', next)
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

module.exports = {
  parseCSV,
  profilerExtraction,
  visualizationProfilerExtraction,
  getResults,
  visualizationCalc,
  visualizationData,
  getSignaturesUser,
  getPublicData,
  upload,
  visualizationDownload,
  visualizationDownloadPublic,
  explorationCalc,
  explorationData,
  submitQueue,
  getQueueResults,
  getVisExample,
  getExposureExample,
  getPublications,
  getImageS3,
  getFileS3,
  to2dArray,
  getRelativePath,
};
