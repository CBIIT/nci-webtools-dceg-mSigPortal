const cluster = require('cluster');
const numCPUs = require('os').cpus().length;
const path = require('path');
const express = require('express');
const config = require('./config.json');
const logger = require('./logger');
const app = express();
const {
  visualizationProfilerExtraction,
  getSummary,
  visualizeR,
  getReferenceSignatureSets,
  getSignatures,
  getPublicDataOptions,
  getPublicData,
  upload,
  getPublicSVG,
  download,
  exploringR,
  getReferenceSignatureData,
  submitQueue,
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

  const server = app.listen(config.port, () => {
    logger.info(
      `msigconfig.portal server running on config.port: ${config.port}`
    );
  });

  server.keepAliveTimeout = 61 * 1000;
  server.headersTimeout = 62 * 1000;
}

app.use(express.static(path.resolve('www')));
app.use('/results', express.static(config.tmppath));
app.use(express.json());

app.use((error, req, res, next) => {
  logger.error(err);
  if (!error.statusCode) error.statusCode = 500;
  return res.status(error.statusCode).json({ error: error.toString() });
});

app.get('/ping', (req, res) => res.send(true));

app.post('/profilerExtraction', visualizationProfilerExtraction);

app.post('/getSummary', getSummary);

app.post('/visualizeR', visualizeR);

app.post('/visualizeR/getReferenceSignatureSets', getReferenceSignatureSets);

app.post('/visualizeR/getSignatures', getSignatures);

app.post('/getPublicDataOptions', getPublicDataOptions);

app.post('/getPublicData', getPublicData);

app.post('/upload', upload);

app.post('/getPublicSVG', getPublicSVG);

app.get('/visualize/download', download);

app.post('/exploringR', exploringR);

app.post('/getReferenceSignatureData', getReferenceSignatureData);

app.post('/queue', submitQueue);
