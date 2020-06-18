const express = require('express');
const app = express();
const cors = require('cors');
var path = require('path');
const logger = require('./logger');
const { port } = require('./config.json');
const { spawn } = require('child_process');
const r = require('r-wrapper')


app.use(cors());
app.use(express.json());
app.use(express.static('www'));

app.get('/ping', (req, res) => res.send(true));

app.get('/api/python', (req, res) => {

  const wrapper = spawn('python3', [
    'api/python/wrapper.py',
    req.query.fn,
    req.query.arg
  ]);

  wrapper.stdout.on('data', (data) => res.send(data.toString()));
  wrapper.stderr.on('data', (data) => console.log(`error: ${data}`));
});

app.post('/r', (req, res) => {

  const data = r(path.resolve(__dirname, 'test.R'),req.body.script,{
    first: req.body.first !== undefined ? req.body.first : 'John',
    last: req.body.last !== undefined ? req.body.last : 'Smith'
  })
  res.send(data.toString())
});


app.listen(port, () => {
  logger.log('info', `Application running on port: ${port}`);
  console.log(`Listening on port ${port}`);
});
