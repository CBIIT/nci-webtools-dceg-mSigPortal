const express = require('express');
const fs = require('fs-extra');
const knex = require('knex');
const logger = require('./services/logger');
const config = require('./config.json');
const serviceRouter = require('./services/extraction/router');

const app = createApp(config);
const server = app.listen(config.extraction.port, () => {
  app.locals.logger.info(
    `mSigPortal extraction service running on config.port: ${config.extraction.port}`
  );
});

function createApp(config) {
  const app = express();

  // register services as app locals
  if (config.database) {
    app.locals.connection = knex({
      client: 'postgres',
      connection: config.database,
    });
  }
  app.locals.logger = logger;
  app.use(serviceRouter);

  // create required folders
  for (let folder of [config.logs.folder, config.results.folder]) {
    fs.mkdirSync(folder, { recursive: true });
  }

  return app;
}
