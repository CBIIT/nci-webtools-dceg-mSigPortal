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
  visualizationWrapper,
  getSignaturesUser,
  upload,
  visualizationDownload,
  visualizationDownloadPublic,
  explorationWrapper,
  submitQueue,
  getQueueResults,
  getVisExample,
  getExposureExample,
  getPublications,
  getImageS3Batch,
  getImageS3,
  getFileS3,
  downloadSession,
  associationWrapper,
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

apiRouter.use(
  '/results',
  express.static(config.results.folder, {
    // setHeaders: (res, path, stat) => {
    //   res.set('Cache-Control', 'max-age=0, must-revalidate');
    // },
  })
);

apiRouter.use(express.json({ limit: '50mb' }));

apiRouter.get('/ping', (req, res) => res.send(true));

apiRouter.post('/profilerExtraction', visualizationProfilerExtraction);

apiRouter.post('/getResults', getResults);

apiRouter.post('/visualizationWrapper', visualizationWrapper);

apiRouter.post('/getSignaturesUser', getSignaturesUser);

apiRouter.post('/upload', upload);

apiRouter.get('/visualization/download', visualizationDownload);

apiRouter.post('/visualization/downloadPublic', visualizationDownloadPublic);

apiRouter.post('/explorationWrapper', explorationWrapper);

apiRouter.post('/queue', submitQueue);

apiRouter.get('/getQueueResults/:id', getQueueResults);

apiRouter.get('/getVisExample/:example', getVisExample);

apiRouter.get('/getExposureExample/:example', getExposureExample);

apiRouter.get('/getPublications', getPublications);

apiRouter.post('/getImageS3Batch', getImageS3Batch);

apiRouter.post('/getImageS3', getImageS3);

apiRouter.post('/getFileS3', getFileS3);

apiRouter.post('/downloadSession', downloadSession);

apiRouter.post('/associationWrapper', associationWrapper);

apiRouter.use((err, req, res, next) => {
  logger.debug(err);
  const { name, message, stack } = err;
  logger.error(err);
  res.status(500).json(`${name}: ${message}`);
});
