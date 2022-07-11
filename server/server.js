const path = require('path');
const express = require('express');
const fs = require('fs');
const knex = require('knex');
const logger = require('./services/logger');
const { apiRouter } = require('./services/router');
const config = require('./config.json');

const app = createApp(config);
const server = app.listen(config.server.port, () => {
  logger.info(
    `msigconfig.portal server running on config.port: ${config.server.port}`
  );
});
server.keepAliveTimeout = 61 * 1000;
server.headersTimeout = 62 * 1000;

function createApp(config) {
  const app = express();

  if (config.database) {
    app.locals.connection = knex({ client: 'postgres', connection: config.database });
  }

  app.use(apiRouter);
  // app.use(express.static(config.server.static));
  // app.use(express.static(path.resolve('www')));

  // create required folders
  for (let folder of [config.logs.folder, config.results.folder]) {
    fs.mkdirSync(folder, { recursive: true });
  }

  return app;
}
