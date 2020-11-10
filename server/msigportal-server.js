const cluster = require('cluster');
const numCPUs = require('os').cpus().length;
const path = require('path');
const express = require('express');
const AWS = require('aws-sdk');
const fs = require('fs');
const config = require('./config.json');
const logger = require('./logger');
const app = express();
const {
  visualizationProfilerExtraction,
  getResultData,
  visualizeR,
  getReferenceSignatureSets,
  getSignatures,
  getPublicDataOptions,
  getPublicData,
  upload,
  download,
  exploringR,
  getReferenceSignatureData,
  submitQueue,
  fetchResults,
  fetchExample,
} = require('./controllers');

if (cluster.isMaster) {
  masterProcess();
} else {
  childProcess();
}

function masterProcess() {
  console.log(`Master ${process.pid} is running`);

  for (let i = 0; i < numCPUs; i++) {
    console.log(`Forking process number ${i}...`);
    cluster.fork();
  }
}

function childProcess() {
  console.log(`Worker ${process.pid} started and finished`);

  const server = app.listen(config.server.port, () => {
    // update aws configuration if supplied
    if (config.aws) {
      AWS.config.update(config.aws);
    }

    // create required folders
    for (let folder of [config.logs.folder, config.results.folder]) {
      fs.mkdirSync(folder, { recursive: true });
    }

    logger.info(
      `msigconfig.portal server running on config.port: ${config.server.port}`
    );
  });

  server.keepAliveTimeout = 61 * 1000;
  server.headersTimeout = 62 * 1000;
}

app.use(express.static(path.resolve('www')));
app.use('/results', express.static(config.results.folder));
app.use('/public', express.static(config.data.database));
app.use('/examples', express.static(path.join(config.data.folder, 'Examples')));
app.use(express.json());

app.use((error, req, res, next) => {
  logger.error(error);
  if (!error.statusCode) error.statusCode = 500;
  res.status(error.statusCode).json(error.message);
});

app.get('/ping', (req, res) => res.send(true));

app.post('/profilerExtraction', visualizationProfilerExtraction);

app.post('/getResultData', getResultData);

app.post('/visualizeR', visualizeR);

app.post('/visualizeR/getReferenceSignatureSets', getReferenceSignatureSets);

app.post('/visualizeR/getSignatures', getSignatures);

app.post('/getPublicDataOptions', getPublicDataOptions);

app.post('/getPublicData', getPublicData);

app.post('/upload', upload);

app.get('/visualize/download', download);

app.post('/exploringR', exploringR);

app.post('/getReferenceSignatureData', getReferenceSignatureData);

app.post('/queue', submitQueue);

app.get('/fetchResults/:id', fetchResults);

app.get('/fetchExample/:folder', fetchExample);
