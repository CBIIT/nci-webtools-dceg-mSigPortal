const cluster = require('cluster');
const numCPUs = require('os').cpus().length;
const path = require('path');
const express = require('express');
const { port } = require('./config.json');
const logger = require('./logger');
const app = express();
const {
  profilerExtraction,
  getSummary,
  visualizeR,
  getReferenceSignatureSets,
  getSignatures,
  getPublicDataOptions,
  getPublicData,
  upload,
  downloadPlotData,
  getSVG,
  getPublicSVG,
  download,
  exploringR,
  getReferenceSignatureData,
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

  const server = app.listen(port, () => {
    logger.info(`msigportal server running on port: ${port}`);
    console.log(`Listening on port ${port}`);
  });

  server.keepAliveTimeout = 61 * 1000;
  server.headersTimeout = 62 * 1000;
}

app.use(express.static(path.resolve('www')));
app.use(express.json());

app.use((err, req, res, next) => {
  logger.info('Unhandled Error:\n' + err.message);
  logger.error(err);
  if (!err.statusCode) err.statusCode = 500;
  res.status(err.statusCode).send(err.message);
});

app.get('/ping', (req, res) => res.send(true));

app.post('/profilerExtraction', profilerExtraction);

app.post('/getSummary', getSummary);

app.post('/visualizeR', visualizeR);

app.post('/visualizeR/getReferenceSignatureSets', getReferenceSignatureSets);

app.post('/visualizeR/getSignatures', getSignatures);

app.post('/getPublicDataOptions', getPublicDataOptions);

app.post('/getPublicData', getPublicData);

app.post('/upload', upload);

app.post('/downloadPlotData', downloadPlotData);

app.post('/getSVG', getSVG);

app.post('/getPublicSVG', getPublicSVG);

app.get('/visualize/download', download);

app.post('/exploringR', exploringR);

app.post('/getReferenceSignatureData', getReferenceSignatureData);

app.post('/visualize/queue', (req, res) => {
  res.json(req.body);
});
