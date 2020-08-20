const express = require('express');
const path = require('path');
const logger = require('./logger');
const { port, tmppath } = require('./config.json');
const { spawn } = require('child_process');
const formidable = require('formidable');
const fs = require('fs');
const rimraf = require('rimraf');
const { v4: uuidv4 } = require('uuid');
const Papa = require('papaparse');
const r = require('r-wrapper');
const AWS = require('aws-sdk');
const app = express();

app.use(express.static(path.resolve('www')));
app.use(express.json());

app.use((err, req, res, next) => {
  logger.info('Unhandled Error:\n' + err.message);
  logger.error(err);
  if (!err.statusCode) err.statusCode = 500;
  res.status(err.statusCode).send(err.message);
});

app.get('/ping', (req, res) => res.send(true));

/**
 * Calls wrapper function for Python and R
 *
 * @param {string} req.body.fn Wrapper function name
 * @param {string} req.body.args Name of function handler (python or R)
 * @return {string} JSON string of function name and arguments
 */
app.post('/api', (req, res) => {
  let stdout = '';
  let stderr = '';
  const wrapper = spawn('python3', [
    'api/wrapper.py',
    req.body.handler,
    JSON.stringify(req.body.params),
  ]);

  wrapper.stderr.on('data', (data) => (stderr += data.toString()));
  wrapper.stderr.on('close', () => {
    console.log('return', { stdout, stderr });
    res.json({ return: stdout, err: stderr });
  });
});

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

async function getSummary(resultsPath) {
  const summaryPath = path.join(resultsPath, 'Summary.txt');
  const statisticsPath = path.join(resultsPath, 'Statistics.txt');
  const matrixPath = path.join(resultsPath, 'Matrix_List.txt');
  const downloadsPath = path.join(resultsPath, 'output');
  let matrixList = [];
  let statistics = '';
  let downloads = [];
  const summary = await parseCSV(summaryPath);

  if (fs.existsSync(matrixPath)) matrixList = await parseCSV(matrixPath);
  if (fs.existsSync(statisticsPath))
    statistics = fs.readFileSync(statisticsPath, 'utf8');

  if (fs.existsSync(downloadsPath)) {
    downloads = fs
      .readdirSync(downloadsPath)
      .filter((file) => file.endsWith('.tar.gz'))
      .map((file) => file);
  }

  return {
    summary: summary,
    statistics: statistics,
    matrixList: matrixList,
    downloads: downloads,
  };
}

app.post('/api/visualize', (req, res) => {
  req.setTimeout(15 * 60 * 1000);
  res.setTimeout(15 * 60 * 1000, () => {
    res.status(504).send('request timed out');
  });

  logger.info('/api/visualize: Spawning Python Process');
  let reqBody = { ...req.body };
  // update paths
  reqBody.outputDir[1] = path.join(tmppath, reqBody.outputDir[1], 'results');
  const args = Object.values(reqBody);
  const cli = args.reduce((params, arg) => [...params, ...arg]);

  logger.debug('/api/visualize: CLI args\n' + cli);
  let stdout = '';
  let stderr = '';
  const wrapper = spawn('python3', [
    'api/python/mSigPortal_Profiler_Extraction.py',
    ...cli,
  ]);

  wrapper.stdout.on('data', (data) => (stdout += data.toString()));
  wrapper.stderr.on('data', (data) => (stderr += data.toString()));
  wrapper.stderr.on('close', async () => {
    const scriptOut = { stdout: stdout, stderr: stderr };
    const resultsPath = path.join(tmppath, reqBody.projectID[1], 'results');
    // logger.debug('STDOUT\n' + scriptOut.stdout);
    // logger.debug('STDERR\n' + scriptOut.stderr);

    if (fs.existsSync(path.join(resultsPath, 'Summary.txt'))) {
      res.json({ ...scriptOut, ...(await getSummary(resultsPath)) });
    } else {
      logger.info('/api/visualize: An Error Occured While Extracting Profiles');
      res.status(500).json({
        msg:
          'An error occured durring profile extraction. Please review your input parameters and try again.',
        ...scriptOut,
      });
    }
  });
});

// read summary file and return plot mapping
app.post('/visualize/summary', async (req, res) => {
  logger.info('/visualize/summary: Retrieving Summary');
  const resultsPath = path.join(tmppath, req.body.projectID, 'results');

  if (fs.existsSync(path.join(resultsPath, 'Summary.txt'))) {
    res.json(await getSummary(resultsPath));
  } else {
    logger.info('/visualize/summary: Summary file not found');
    res.status(500).json({
      msg: 'Summary file not found',
    });
  }
});

app.post('/api/visualizeR', (req, res) => {
  logger.info('/api/visualizeR: function ' + req.body.fn);
  console.log('args', req.body);
  const savePath = path.join(
    tmppath,
    req.body.projectID,
    'results',
    req.body.fn
  );

  if (!fs.existsSync(savePath)) {
    fs.mkdirSync(savePath);
  }
  try {
    const wrapper = r('api/R/visualizeWrapper.R', req.body.fn, {
      ...req.body.args,
      projectID: req.body.projectID,
      pythonOutput: path.join(tmppath, req.body.projectID, 'results/output'),
      savePath: savePath,
    });

    const { stdout, output } = JSON.parse(wrapper);
    console.log('wrapper return', JSON.parse(wrapper));

    res.json({
      debugR: stdout,
      output: output,
    });
  } catch (err) {
    logger.info('/api/visualizeR: An error occured');
    logger.error(err);
    res.status(500).json(err.message);
  }
});

app.post('/api/visualizeR/getReferenceSignatureSets', (req, res) => {
  logger.info('/api/visualizeR/getReferenceSignatureSets: Calling R Wrapper');
  console.log('args', req.body);

  try {
    const list = r('api/R/visualizeWrapper.R', 'getReferenceSignatureSets', [
      req.body.profileType,
    ]);

    console.log('SignatureReferenceSets', list);

    res.json(list);
  } catch (err) {
    logger.info('/api/visualizeR/getReferenceSignatureSets: An error occured');
    logger.error(err);
    res.status(500).json(err.message);
  }
});

app.post('/api/visualizeR/getSignatures', (req, res) => {
  logger.info('/api/visualizeR/getSignatures: Calling R Wrapper');

  try {
    const list = r('api/R/visualizeWrapper.R', 'getSignatures', [
      req.body.profileType,
      req.body.signatureSetName,
    ]);

    console.log('signatures', list);

    res.json(list);
  } catch (err) {
    logger.info('/api/visualizeR/getSignatures: An error occured');
    logger.error(err);
    res.status(500).json(err.message);
  }
});

app.post('/visualize/upload', (req, res, next) => {
  const projectID = uuidv4();
  const form = formidable({
    uploadDir: path.join(tmppath, projectID),
    multiples: true,
  });

  logger.info(`/VISUALIZE/UPLOAD: Request Project ID:${projectID}`);

  if (!fs.existsSync(form.uploadDir)) {
    fs.mkdirSync(form.uploadDir);
  } else {
    rimraf(form.uploadDir, () => {
      fs.mkdirSync(form.uploadDir);
    });
  }

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
});

app.post('/visualize/txt', (req, res) => {
  const txtPath = req.body.path;
  if (txtPath.indexOf(tmppath) == 0) {
    const s = fs.createReadStream(txtPath);

    s.on('open', () => {
      res.set('Content-Type', 'text/plain');
      s.pipe(res);
      logger.debug(`/visualize/txt: Serving ${txtPath}`);
    });
    s.on('error', () => {
      res.set('Content-Type', 'text/plain');
      res.status(500).end('Not found');
      logger.info(`/visualize/txt: Error retrieving ${txtPath}`);
    });
  } else {
    logger.info('traversal error');
    res.status(500).end('Not found');
  }
});

app.post('/visualize/svg', (req, res) => {
  const svgPath = req.body.path;
  if (svgPath.indexOf(tmppath) == 0) {
    const s = fs.createReadStream(svgPath);

    s.on('open', () => {
      res.set('Content-Type', 'image/svg+xml');
      s.pipe(res);
      logger.debug(`/VISUALIZE/SVG: Serving ${svgPath}`);
    });
    s.on('error', () => {
      res.set('Content-Type', 'text/plain');
      res.status(500).end('Not found');
      logger.info(`/VISUALIZE/SVG: Error retrieving ${svgPath}`);
    });
  } else {
    logger.info('traversal error');
    res.status(500).end('Not found');
  }
});

// download generated result outputs
app.get('/visualize/download', (req, res) => {
  logger.info(`/visualize/download: id:${req.query.id} file:${req.query.file}`);
  const file = path.join(
    tmppath,
    req.query.id,
    'results/output',
    req.query.file
  );
  if (file.indexOf(tmppath) == 0) {
    res.download(file);
  } else {
    logger.info('traversal error');
    res.status(500).end('Not found');
  }
});

app.post('/visualize/queue', (req, res) => {
  res.json(req.body);
});

app.listen(port, () => {
  logger.info(`msigportal server running on port: ${port}`);
  console.log(`Listening on port ${port}`);
});
