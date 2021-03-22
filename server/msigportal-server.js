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
  visualizeR,
  getReferenceSignatureSets,
  getSignaturesR,
  getSignaturesUser,
  getPublicData,
  upload,
  download,
  exploringR,
  getReferenceSignatureData,
  submitQueue,
  getQueueResults,
  getVisExample,
  getSignatureNames,
  getSampleNames,
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

// serve public folder during local development
if (process.env.NODE_ENV !== 'production')
  // app.use(express.static(config.server.static));
  app.use(express.static(path.resolve('www')));

apiRouter.use('/results', express.static(config.results.folder));
apiRouter.use('/public', express.static(config.data.localDatabase));
apiRouter.use(express.json());

apiRouter.use((error, req, res, next) => {
  logger.error(error);
  if (!error.statusCode) error.statusCode = 500;
  res.status(error.statusCode).json(error.message);
});

apiRouter.get('/ping', (req, res) => res.send(true));

apiRouter.post('/profilerExtraction', visualizationProfilerExtraction);

apiRouter.post('/getResults', getResults);

apiRouter.post('/visualizeR', visualizeR);

apiRouter.post('/getReferenceSignatureSets', getReferenceSignatureSets);

apiRouter.post('/getSignaturesR', getSignaturesR);

apiRouter.post('/getSignaturesUser', getSignaturesUser);

apiRouter.post('/getPublicData', getPublicData);

apiRouter.post('/upload', upload);

apiRouter.get('/visualize/download', download);

apiRouter.post('/exploringR', exploringR);

apiRouter.post('/getReferenceSignatureData', getReferenceSignatureData);

apiRouter.post('/queue', submitQueue);

apiRouter.get('/getQueueResults/:id', getQueueResults);

apiRouter.get('/getVisExample/:example', getVisExample);

apiRouter.post('/getSignatureNames', getSignatureNames);

apiRouter.post('/getSampleNames', getSampleNames);

apiRouter.get('/getExposureExample/:example', getExposureExample);

apiRouter.get('/getPublications', getPublications);

apiRouter.post('/getImageS3', getImageS3);

apiRouter.post('/getFileS3', getFileS3);
