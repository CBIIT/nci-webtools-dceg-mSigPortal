const cluster = require('cluster');
const numCPUs = require('os').cpus().length;
const path = require('path');
const express = require('express');
const fs = require('fs');
const config = require('./config.json');
const logger = require('./logger');
const compression = require('compression');

const app = express();
const apiRouter = express.Router();

app.use(compression());
app.use('/api', apiRouter);

const {
  visualizationProfilerExtraction,
  getResults,
  visualizationCalc,
  visualizationData,
  getSignaturesUser,
  getPublicData,
  upload,
  visualizationDownload,
  visualizationDownloadPublic,
  explorationCalc,
  explorationData,
  submitQueue,
  getQueueResults,
  getVisExample,
  getExposureExample,
  getPublications,
  getImageS3,
  getFileS3,
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

// app.use(express.static(config.server.static));
app.use(express.static(path.resolve('www')));

apiRouter.use('/results', express.static(config.results.folder));
// apiRouter.use('/public', express.static(config.data.localData));
apiRouter.use(express.json());

apiRouter.use((error, req, res, next) => {
  const { name, message, stack } = error;
  logger.error({ message, stack });
  response.status(500).json(`${name}: ${message}`);
});

apiRouter.get('/ping', (req, res) => res.send(true));

apiRouter.post('/profilerExtraction', visualizationProfilerExtraction);

apiRouter.post('/getResults', getResults);

apiRouter.post('/visualizationCalc', visualizationCalc);

apiRouter.post('/visualizationData', visualizationData);

apiRouter.post('/getSignaturesUser', getSignaturesUser);

apiRouter.post('/getPublicData', getPublicData);

apiRouter.post('/upload', upload);

apiRouter.get('/visualization/download', visualizationDownload);

apiRouter.post('/visualization/downloadPublic', visualizationDownloadPublic);

apiRouter.post('/explorationCalc', explorationCalc);

apiRouter.post('/explorationData', explorationData);

apiRouter.post('/queue', submitQueue);

apiRouter.get('/getQueueResults/:id', getQueueResults);

apiRouter.get('/getVisExample/:example', getVisExample);

apiRouter.get('/getExposureExample/:example', getExposureExample);

apiRouter.get('/getPublications', getPublications);

apiRouter.post('/getImageS3', getImageS3);

apiRouter.post('/getFileS3', getFileS3);
