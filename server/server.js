const cluster = require('cluster');
const numCPUs = require('os').cpus().length;
const path = require('path');
const express = require('express');
const fs = require('fs');
const knex = require('knex');
const config = require('./config.json');
const logger = require('./services/logger');
const { apiRouter } = require('./services/api');

const app = express();

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

  if (config.database) {
    app.locals.connection = knex({ client: 'postgres', connection: config.database });
  }

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
// app.use(express.static(path.resolve('www')));

app.use(apiRouter);
