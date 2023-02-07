const path = require('path');
const express = require('express');
const fs = require('fs');
const knex = require('knex');
const logger = require('./services/logger');
const { createDatabaseCache } = require('./services/cache');
const { apiRouter } = require('./routes/router');
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
    app.locals.connection = knex({
      client: 'postgres',
      connection: config.database,
    });
  }
  app.locals.sqlite = (userId, db) =>
    knex({
      client: 'better-sqlite3',
      connection: () => ({
        filename: path.join(config.results.folder, userId, `${db}.sqlite3`),
      }),
      useNullAsDefault: true,
    });
  app.locals.cache = createDatabaseCache(app.locals.connection, 'cache');
  app.locals.cache.initialize();

  app.use(apiRouter);
  // app.use(express.static(config.server.static));
  // app.use(express.static(path.resolve('www')));

  // create required folders
  for (let folder of [config.logs.folder, config.results.folder]) {
    fs.mkdirSync(folder, { recursive: true });
  }

  return app;
}
