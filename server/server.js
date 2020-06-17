const express = require('express');
const cors = require('cors');
const path = require('path');
const logger = require('./logger');
const { port } = require('./config.json');
const { spawn } = require('child_process');

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
 * @param {string} req.body.args JSON string of function params
 * @return {string} JSON string of function output
 */
app.post('/api', (req, res) => {
  let stdout = '';
  let stderr = '';
  const wrapper = spawn('python3', [
    'api/wrapper.py',
    req.body.fn,
    req.body.args,
  ]);

  wrapper.stdout.on('data', (data) => (stdout += data.toString()));
  wrapper.stderr.on('data', (data) => (stderr += data.toString()));
  wrapper.stderr.on('close', () => {
    console.log('return', { stdout, stderr });
    res.send({ return: stdout, err: stderr });
  });
});

app.listen(port, () => {
  logger.log('info', `Application running on port: ${port}`);
  console.log(`Listening on port ${port}`);
});
