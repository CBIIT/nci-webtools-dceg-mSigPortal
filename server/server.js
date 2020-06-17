const express = require('express');
const app = express();
const cors = require('cors');
var path = require('path');
const logger = require('./logger');
const { port } = require('./config.json');
const { spawn } = require('child_process');

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

app.get('/api/r', (req, res) => {
  
   //Validate args and check for optional/missing args here. Assign default for optional args

  const wrapper= spawn('python3', [
    "api/python/r-wrapper.py",
    req.query.fn,
    req.query.first,
    req.query.second 
  ]);
 
  wrapper.stdout.on('data', (data) => res.send(data.toString()));
  wrapper.stderr.on('data', (data) => console.log(`error: ${data}`));
});


app.listen(port, () => {
  logger.log('info', `Application running on port: ${port}`);
  console.log(`Listening on port ${port}`);
});
