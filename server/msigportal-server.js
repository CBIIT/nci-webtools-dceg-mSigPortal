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
var cron = require('node-cron');

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

  wrapper.stdout.on('data', (data) => (stdout += data.toString()));
  wrapper.stderr.on('data', (data) => (stderr += data.toString()));
  wrapper.stderr.on('close', () => {
    console.log('return', { stdout, stderr });
    res.json({ return: stdout, err: stderr });
  });
});

app.post('/api/visualize', (req, res) => {
  logger.info('/API/VISUALIZE: Spawning Python Process');
  let reqBody = { ...req.body };
  // resolve paths
  const root = 'server/tmp';
  const resultsDir = path.join(reqBody.outputDir[1], 'results');
  reqBody.inputFile[1] = path.resolve(root, reqBody.inputFile[1]);
  reqBody.outputDir[1] = path.resolve(root, resultsDir);
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
    console.log('return', { stdout, stderr });
    if (stderr.length) {
      res.json({ return: stdout, err: stderr });
    } else {
      let plotDir = '';
      if (reqBody.inputFormat[1].includes('catalog'))
        plotDir = path.join(resultsDir, 'svg');
      else plotDir = path.join(resultsDir, 'output', 'plots', 'svg');
      const plots = fs.readdirSync(path.resolve(root, plotDir));
      console.log(plots);

      res.json({ plotDir: plotDir, plots: plots, return: stdout, err: stderr });
    }
  });
});

app.post('/upload', (req, res, next) => {
  const projectID = uuidv4();
  const form = formidable({ uploadDir: path.resolve(tmppath, projectID) });

  logger.info(`/UPLOAD: Request UUID:${projectID}`);

  if (!fs.existsSync(form.uploadDir)) {
    fs.mkdirSync(form.uploadDir);
  } else {
    rimraf(form.uploadDir, () => {
      fs.mkdirSync(form.uploadDir);
    });
  }

  app.post('/svg', (req, res) => {
    const root = 'server/tmp';
    const svgPath = path.resolve(root, req.body.path);
    console.log(svgPath);
    var s = fs.createReadStream(svgPath);
    s.on('open', () => {
      res.set('Content-Type', 'image/svg+xml');
      console.log('open');
      s.pipe(res);
    });
    s.on('error', () => {
      res.set('Content-Type', 'text/plain');
      console.log('error');
      res.status(404).end('Not found');
    });
  });

  form.parse(req);
  form.on('file', (field, file) => {
    const uploadPath = path.join(form.uploadDir, file.name);
    fs.rename(file.path, uploadPath, (err) => {
      if (err) {
        logger.error(
          `/UPLOAD: Failed to upload file: ${file.name}` + '\n' + err
        );
        res.status(404).json({
          msg: `/UPLOAD: Failed to upload file: ${file.name}`,
          error: err,
        });
      } else {
        logger.info(`/UPLOAD: Successfully uploaded file: ${file.name}`);
        res.json({
          projectID: projectID,
          filePath: path.join(projectID, file.name),
          // filePath: uploadPath,
          // outputDir: form.uploadDir,
        });
      }
    });
  });
  form.on('error', function (err) {
    logger.error('/UPLOAD: An error occured\n' + err);
    res.status(500).json({
      msg: '/UPLOAD: An error has occured',
      err: err,
    });
  });
});

cron.schedule('50 7 * * *', function () {
  const process = spawn(
    'find local/content/analysistools/public_html/apps/msigportal/tmp -mindepth 1 -mtime +14 -exec rm {} ; >>var/log/msigportal-cron.log 2>&1',
    [],
    { shell: true }
  );
  process.stderr.on('data', (data) => console.log(data.toString()));
});

cron.schedule('45 7 * * *', function () {
  const process = spawn(
    'find  local/content/analysistools/public_html/apps/msigportal/logs -mindepth 1 -mtime +60 -exec rm {} ; >>var/log/msigportal-cron.log 2>&1',
    [],
    { shell: true }
  );
  process.stderr.on('data', (data) => console.log(data.toString()));
});

app.listen(port, () => {
  logger.info(`msigportal server running on port: ${port}`);
  console.log(`Listening on port ${port}`);
});
