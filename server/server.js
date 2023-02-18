import path from 'path';
import express from 'express';
import fs from 'fs';
import knex from 'knex';
import logger from './services/logger.js';
import { createDatabaseCache } from './services/cache.js';
import { apiRouter } from './routes/router.js';
import config from './config.json' assert { type: 'json' };

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
      connection: {
        filename: path.join(config.folders.output, userId, `${db}.sqlite3`),
      },
      useNullAsDefault: true,
    });
  app.locals.cache = createDatabaseCache(app.locals.connection, 'cache');
  app.locals.cache.initialize();

  app.use(apiRouter);

  // create required folders
  for (let folder of [config.logs.folder, config.results.folder]) {
    fs.mkdirSync(folder, { recursive: true });
  }

  return app;
}
