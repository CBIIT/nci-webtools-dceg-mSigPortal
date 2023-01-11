const express = require('express');
const fs = require('fs-extra');
const knex = require('knex');
const logger = require('./services/logger');
const config = require('./config.json');
const apiRouter = require('./services/extraction/router');

const app = createApp(config);
const server = app.listen(config.extraction.port, () => {
  logger.info(
    `mSigPortal extraction service running on config.port: ${config.extraction.port}`
  );
});
server.keepAliveTimeout = 61 * 1000;
server.headersTimeout = 62 * 1000;

function createApp(config) {
  const app = express();

  if (config.database) {
    app.locals.connection = knex({
      client: 'postgres',
      connection: config.database,
    });
  }

  app.use(apiRouter);

  // create required folders
  for (let folder of [config.logs.folder, config.results.folder]) {
    fs.mkdirSync(folder, { recursive: true });
  }

  return app;
}
