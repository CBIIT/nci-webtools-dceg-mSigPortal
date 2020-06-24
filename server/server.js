const express = require('express');
const cors = require('cors');
const path = require('path');
const logger = require('./logger');
const { port } = require('./config.json');
const { spawn } = require('child_process');
const r = require('r-wrapper')
const cron = require("node-cron");

const app = express();

app.use(cors());
app.use(express.json());
app.use(express.static('www'));
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
    res.send({ return: stdout, err: stderr });
  });
});

cron.schedule("50 7 * * *", function() {
  const process = spawn('find local/content/analysistools/public_html/apps/msigportal/tmp -mindepth 1 -mtime +14 -exec rm {} \; >>var/log/msigportal-cron.log 2>&1', [],{shell:true})
  process.stderr.on('data',(data) => console.log(data.toString()))
});

cron.schedule("45 7 * * *", function() {
  const process = spawn('find  local/content/analysistools/public_html/apps/msigportal/logs -mindepth 1 -mtime +60 -exec rm {} \; >>var/log/msigportal-cron.log 2>&1', [],{shell:true})
  process.stderr.on('data',(data) => console.log(data.toString()))
});

app.listen(port, () => {
  logger.log('info', `Application running on port: ${port}`);
  console.log(`Listening on port ${port}`);
});
