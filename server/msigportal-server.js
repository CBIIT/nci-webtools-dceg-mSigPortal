const cluster = require('cluster');
const numCPUs = require('os').cpus().length;
const path = require('path');
const express = require('express');
const AWS = require('aws-sdk');
const fs = require('fs');
const config = require('./config.json');
const logger = require('./logger');

const app = express();
const apiRouter = express.Router();
app.use('/api', apiRouter);

const {
  visualizationProfilerExtraction,
  getResultData,
  visualizeR,
  getReferenceSignatureSets,
  getSignaturesR,
  getSignaturesUser,
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

// serve public folder during local development
if (process.env.NODE_ENV !== 'production')
  // app.use(express.static(config.server.static));
  app.use(express.static(path.resolve('www')));

apiRouter.use('/results', express.static(config.results.folder));
apiRouter.use('/public', express.static(config.data.database));
apiRouter.use(express.json());

apiRouter.use((error, req, res, next) => {
  logger.error(error);
  if (!error.statusCode) error.statusCode = 500;
  res.status(error.statusCode).json(error.message);
});

apiRouter.get('/ping', (req, res) => res.send(true));

apiRouter.post('/profilerExtraction', visualizationProfilerExtraction);

apiRouter.post('/getResultData', getResultData);

apiRouter.post('/visualizeR', visualizeR);

apiRouter.post('/getReferenceSignatureSets', getReferenceSignatureSets);

apiRouter.post('/getSignaturesR', getSignaturesR);

apiRouter.post('/getSignaturesUser', getSignaturesUser);

apiRouter.post('/getPublicDataOptions', getPublicDataOptions);

apiRouter.post('/getPublicData', getPublicData);

apiRouter.post('/upload', upload);

apiRouter.get('/visualize/download', download);

apiRouter.post('/exploringR', exploringR);

apiRouter.post('/getReferenceSignatureData', getReferenceSignatureData);

apiRouter.post('/queue', submitQueue);

apiRouter.get('/fetchResults/:id', fetchResults);

apiRouter.get('/fetchExample/:example', fetchExample);
