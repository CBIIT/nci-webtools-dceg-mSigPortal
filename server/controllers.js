const path = require('path');
const logger = require('./logger');
const { tmppath, datapath } = require('./config.json');
const { spawn } = require('child_process');
const formidable = require('formidable');
const fs = require('fs');
const { v4: uuidv4 } = require('uuid');
const Papa = require('papaparse');
const r = require('r-wrapper').async;
const AWS = require('aws-sdk');

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
    newPaths[key] = fullPath.replace(path.resolve(tmppath) + '/', '');
  });
  return newPaths;
}

async function profilerExtraction(req, res, next) {
  req.setTimeout(15 * 60 * 1000);
  res.setTimeout(15 * 60 * 1000, () => {
    res.status(504).send('request timed out');
  });

  logger.info('/profilerExtraction: Spawning Python Process');
  let reqBody = { ...req.body };
  // update paths
  reqBody.outputDir[1] = path.join(tmppath, reqBody.outputDir[1], 'results');
  const args = Object.values(reqBody);
  const cli = args.reduce((params, arg) => [...params, ...arg]);

  logger.debug('/profilerExtraction: CLI args\n' + cli);
  let stdout = '';
  let stderr = '';
  const wrapper = spawn('python3', [
    'services/python/mSigPortal_Profiler_Extraction.py',
    ...cli,
  ]);

  wrapper.stdout.on('data', (data) => (stdout += data.toString()));
  wrapper.stderr.on('data', (data) => (stderr += data.toString()));
  wrapper.stderr.on('close', async () => {
    const scriptOut = { stdout: stdout, stderr: stderr };
    const resultsPath = path.join(tmppath, reqBody.projectID[1], 'results');
    // logger.debug('STDOUT\n' + scriptOut.stdout);
    // logger.debug('STDERR\n' + scriptOut.stderr);

    if (fs.existsSync(path.join(resultsPath, 'svg_files_list.txt'))) {
      res.json({ ...scriptOut, ...(await getSummaryFiles(resultsPath)) });
    } else {
      logger.info(
        '/profilerExtraction: An Error Occured While Extracting Profiles'
      );
      res.status(500).json({
        msg:
          'An error occured durring profile extraction. Please review your input parameters and try again.',
        ...scriptOut,
      });
    }
  });
}

async function getSummary(req, res, next) {
  logger.info('/getSummary: Retrieving Summary');
  console.log('summary', req.body);
  const resultsPath = path.join(tmppath, req.body.projectID, 'results');

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
    tmppath,
    req.body.projectID,
    'results',
    req.body.fn
  );

  fs.mkdirSync(savePath, { recursive: true });

  try {
    const wrapper = await r('services/R/visualizeWrapper.R', req.body.fn, {
      ...req.body.args,
      projectID: req.body.projectID,
      pythonOutput: path.join(tmppath, req.body.projectID, 'results/output'),
      savePath: savePath,
      dataPath: path.join(datapath, 'signature_visualization/'),
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
  logger.info('/visualizeR/getReferenceSignatureSets: Calling R Wrapper');
  console.log('args', req.body);

  try {
    const list = await r(
      'services/R/visualizeWrapper.R',
      'getReferenceSignatureSets',
      [req.body.profileType, path.join(datapath, 'signature_visualization/')]
    );

    // console.log('SignatureReferenceSets', list);

    res.json(list);
  } catch (err) {
    logger.info('/visualizeR/getReferenceSignatureSets: An error occured');
    logger.error(err);
    res.status(500).json(err.message);
  }
}

async function getSignatures(req, res, next) {
  logger.info('/visualizeR/getSignatures: Calling R Wrapper');

  try {
    const list = await r('services/R/visualizeWrapper.R', 'getSignatures', [
      req.body.profileType,
      req.body.signatureSetName,
      path.join(datapath, 'signature_visualization/'),
    ]);

    // console.log('signatures', list);

    res.json(list);
  } catch (err) {
    logger.info('/visualizeR/getSignatures: An error occured');
    logger.error(err);
    res.status(500).json(err.message);
  }
}

async function getPublicDataOptions(req, res, next) {
  logger.info('/visualize/getPublicOptions: Request');
  try {
    const list = await r(
      'services/R/visualizeWrapper.R',
      'getPublicDataOptions',
      [path.join(datapath, 'signature_visualization/')]
    );

    res.json(JSON.parse(list));
    logger.info('/visualize/getPublicOptions: Success');
  } catch (err) {
    logger.info('/visualize/getPublicOptions: An error occured');
    logger.error(err);
    res.status(500).json(err.message);
  }
}

async function getPublicData(req, res, next) {
  logger.info('/visualize/getPublicOptions: Calling R Wrapper');
  try {
    const list = await r('services/R/visualizeWrapper.R', 'getPublicData', [
      req.body.study,
      req.body.cancerType,
      req.body.experimentalStrategy,
      path.join(datapath, 'signature_visualization/'),
    ]);
    logger.info('/visualize/getPublicOptions: Complete');

    const projectID = uuidv4();
    const saveDir = path.join(tmppath, projectID);

    fs.mkdirSync(saveDir);

    res.json({ svgList: JSON.parse(list), projectID: projectID });
  } catch (err) {
    logger.info('/visualize/getPublicOptions: An error occured');
    logger.error(err);
    res.status(500).json(err.message);
  }
}

function upload(req, res, next) {
  const projectID = uuidv4();
  const form = formidable({
    uploadDir: path.join(tmppath, projectID),
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

function getSVG(req, res, next) {
  const svgPath = path.resolve(req.body.path);
  if (svgPath.indexOf(path.resolve(tmppath)) == 0) {
    const s = fs.createReadStream(svgPath);

    s.on('open', () => {
      res.set('Content-Type', 'image/svg+xml');
      s.pipe(res);
      logger.debug(`/getSVG: Serving ${svgPath}`);
    });
    s.on('error', () => {
      res.set('Content-Type', 'text/plain');
      res.status(500).end('Not found');
      logger.info(`/getSVG: Error retrieving ${svgPath}`);
    });
  } else {
    logger.info('traversal error');
    res.status(500).end('Not found');
  }
}

function getPublicSVG(req, res, next) {
  const svgPath = path.resolve(req.body.path);
  if (svgPath.indexOf(path.resolve(datapath)) == 0) {
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
    path.join(tmppath, req.query.id, 'results/output', req.query.file)
  );
  if (file.indexOf(path.resolve(tmppath)) == 0) {
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
  const savePath = path.join(tmppath, projectID, 'results', req.body.fn);

  fs.mkdirSync(savePath, { recursive: true });

  const wrapper = await r('services/R/exploringWrapper.R', req.body.fn, {
    ...req.body.args,
    projectID: projectID,
    pythonOutput: path.join(tmppath, req.body.projectID, 'results/output'),
    savePath: savePath,
    dataPath: path.join(datapath, 'signature_visualization/'),
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
      dataPath: path.join(datapath, 'signature_visualization/'),
    }
  ).catch(next);

  logger.info('/getReferenceSignatures: Success');
  return res.json(JSON.parse(data));
}

async function submitQueue(res, req, next) {}

module.exports = {
  profilerExtraction,
  getSummary,
  visualizeR,
  getReferenceSignatureSets,
  getSignatures,
  getPublicDataOptions,
  getPublicData,
  upload,
  getSVG,
  getPublicSVG,
  download,
  exploringR,
  getReferenceSignatureData,
  submitQueue,
};
