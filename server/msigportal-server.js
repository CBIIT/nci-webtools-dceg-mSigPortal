const express = require('express');
const cors = require('cors');
const path = require('path');
const logger = require('./logger');
const { port, tmppath } = require('./config.json');
const { spawn } = require('child_process');
const formidable = require('formidable');
const fs = require('fs');
const rimraf = require('rimraf');
const { v4: uuidv4 } = require('uuid');
const Papa = require('papaparse');
const tar = require('tar');
const AWS = require('aws-sdk');

const app = express();

app.use(cors());
app.use(express.static(path.resolve('www')));
app.use(express.json());

app.use((err, req, res, next) => {
  console.error(err.message);
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

app.post('/api/visualize', (req, res) => {
  req.setTimeout(600000);
  logger.info('/API/VISUALIZE: Spawning Python Process');
  let reqBody = { ...req.body };
  // update paths
  const resultsDir = path.join(reqBody.outputDir[1], 'results');
  reqBody.inputFile[1] = path.join(tmppath, reqBody.inputFile[1]);
  reqBody.outputDir[1] = path.join(tmppath, resultsDir);
  const args = Object.values(reqBody);
  const cli = args.reduce((params, arg) => [...params, ...arg]);

  logger.debug('/API/VISUALIZE: CLI args\n' + cli);
  let stdout = '';
  let stderr = '';
  const wrapper = spawn('python3', [
    'api/python/mSigPortal_Profiler_Extraction.py',
    ...cli,
  ]);

  wrapper.stdout.on('data', (data) => (stdout += data.toString()));
  wrapper.stderr.on('data', (data) => (stderr += data.toString()));
  wrapper.stderr.on('close', () => {
    const scriptOut = { stdout: stdout, stderr: stderr };

    logger.debug(`STDOUT ${scriptOut.stdout}`);
    logger.debug(`STDERR: ${scriptOut.stderr}`);

    if (
      fs.existsSync(
        path.join(tmppath, reqBody.projectID[1], 'results', 'Summary.txt')
      )
    ) {
      logger.info('/API/VISUALIZE: Profile Extraction Succeeded');
      res.json({ ...scriptOut });
    } else {
      logger.error(
        '/API/VISUALIZE: An Error Occured While Extracting Profiles'
      );
      res.status(500).json({
        msg:
          'An error occured durring profile extraction. Please review your input parameters and try again.',
        ...scriptOut,
      });
    }
  });
});

app.post('/visualize/upload', (req, res, next) => {
  const projectID = uuidv4();
  const form = formidable({ uploadDir: path.join(tmppath, projectID) });

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
    fs.rename(file.path, uploadPath, (err) => {
      if (err) {
        logger.error(
          `/UPLOAD: Failed to upload file: ${file.name}` + '\n' + err
        );
        res.status(404).json({
          msg: `Failed to upload file: ${file.name}<br><br>
        Review your selected File Type and try again.`,
          error: err,
        });
      } else {
        logger.info(`/UPLOAD: Successfully uploaded file: ${file.name}`);
        res.json({
          projectID: projectID,
          filePath: path.join(projectID, file.name),
        });
      }
    });
  });
  form.on('error', function (err) {
    logger.error('/UPLOAD: An error occured\n' + err);
    res.status(500).json({
      msg: 'An error occured while trying to upload',
      err: err,
    });
  });
});

app.post('/visualize/svg', (req, res) => {
  const svgPath = req.body.path;
  const s = fs.createReadStream(svgPath);

  s.on('open', () => {
    res.set('Content-Type', 'image/svg+xml');
    s.pipe(res);
    logger.debug(`/VISUALIZE/SVG: Serving ${svgPath}`);
  });
  s.on('error', () => {
    res.set('Content-Type', 'text/plain');
    res.status(500).end('Not found');
    logger.error(`/VISUALIZE/SVG: Error retrieving ${svgPath}`);
  });
});

// read summary file and return plot mapping
app.post('/visualize/summary', (req, res) => {
  logger.info('/VISUALIZE/SUMMARY: Reading Results Summary');
  const projectID = req.body.projectID;
  const summaryPath = path.join(tmppath, projectID, 'results', 'Summary.txt');

  fs.readFile(summaryPath, 'utf8', (err, data) => {
    if (err) {
      logger.error('/VISUALIZE/SUMMARY: Error Reading Results Summary');
      logger.error(err);
      res.status(500).json(err);
    }
    Papa.parse(data.trim(), {
      header: true,
      complete: (parsed) => {
        res.json(parsed.data);
      },
    });
  });
});

app.post('/visualize/results', (req, res) => {
  // add error and path traversal checks
  logger.info('/VISUALIZE/RESULTS: Creating Results Archive');
  const projectID = req.body.projectID;
  const archivePath = path.join(tmppath, projectID, 'results.tgz');

  tar
    .create({ sync: true, gzip: true, cwd: path.join(tmppath, projectID) }, [
      'results',
    ])
    .pipe(fs.createWriteStream(archivePath));

  res.json({ projectID: projectID });
});

app.get('/visualize/download', (req, res) => {
  // add logging and checks for path traversal
  const archivePath = path.join(tmppath, req.query.id, 'results.tgz');
  res.download(archivePath, 'results.tgz');
});

app.post('/visualize/queue', (req, res) => {
  res.json(req.body);
});

app.listen(port, () => {
  logger.info(`msigportal server running on port: ${port}`);
  console.log(`Listening on port ${port}`);
});
