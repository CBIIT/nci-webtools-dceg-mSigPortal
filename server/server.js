import path from 'path';
import express from 'express';
import knex from 'knex';
import { validateEnvironment } from './services/environment.js';
import { createLogger } from './services/logger.js';
import { createWebApi, createPublicApi } from './routes/router.js';
import { isMainModule } from './services/utils.js';
import { createDatabaseCache } from './services/cache.js';

// if this module is the main module, start the app
if (isMainModule(import.meta)) {
  const env = process.env;
  validateEnvironment(env);
  main(env);
}

/**
 * Creates and starts an express app given a specified environment.
 * @param {object} env
 * @returns {http.Server} a node http server
 */
export function main(env) {
  const { APP_PORT, APP_NAME, SERVER_TIMEOUT } = env;
  const serverTimeout = +SERVER_TIMEOUT || 1000 * 60;
  const app = createApp(env);
  const server = app.listen(APP_PORT, () => {
    app.locals.logger.info(`${APP_NAME} started on port ${APP_PORT}`);
  });
  server.setTimeout(serverTimeout);
  return server;
}

function createApp(env) {
  const { APP_NAME, LOG_LEVEL } = env;
  const app = express();

  // if behind a proxy, use the first x-forwarded-for address as the client's ip address
  app.set('trust proxy', true);
  app.set('json spaces', 2);
  app.set('x-powered-by', false);

  // register services as app locals
  app.locals.logger = createLogger(APP_NAME, LOG_LEVEL);
  app.locals.connection = knex({
    client: 'postgres',
    connection: {
      host: env.POSTGRES_HOST,
      port: env.POSTGRES_PORT,
      user: env.POSTGRES_USER,
      password: env.POSTGRES_PASS,
      database: env.POSTGRES_DB,
    },
  });
  app.locals.sqlite = (id, filename) =>
    knex({
      client: 'better-sqlite3',
      connection: {
        filename: path.join(env.OUTPUT_FOLDER, id, `${filename}.sqlite3`),
      },
      useNullAsDefault: true,
    });
  app.locals.cache = createDatabaseCache(app.locals.connection, 'cache');
  app.locals.cache.initialize();

  app.use('/web', createWebApi(env));
  app.use('/api', createPublicApi(env));

  return app;
}
