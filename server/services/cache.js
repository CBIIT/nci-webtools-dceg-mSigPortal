export function createDatabaseCache(connection, tableName = 'cache') {
  return {
    initialize: async () => {
      if (!(await connection.schema.hasTable(tableName))) {
        await connection.schema.createTable(tableName, (table) => {
          table.text('key').unique();
          table.json('value');
        });
      }
    },
    get: async (key) => {
      const result = await connection(tableName).where({ key }).first();
      return result ? result.value : null;
    },
    set: async (key, value) => {
      await connection(tableName)
        .insert({ key, value })
        .onConflict('key')
        .merge();
    },
    clear: async (key) => {
      await connection(tableName).where({ key }).del();
    },
    reset: async () => {
      await connection(tableName).truncate();
    },
  };
}

// in-flight responses by cache key; concurrent identical requests in this process share one computation
const pending = new Map();

export function createCacheMiddleware(getCacheKey) {
  return async (req, res, next) => {
    try {
      const key = await getCacheKey(req);
      const cache = req.app.locals.cache;
      const value = await cache.get(key);
      if (value) {
        return res.json(value);
      }

      // undefined means the in-flight request failed; retry with one new leader at a time
      while (pending.has(key)) {
        const body = await pending.get(key);
        if (body !== undefined) return res.json(body);
      }

      let settle;
      const promise = new Promise((resolve) => (settle = resolve));
      pending.set(key, promise);
      const finish = (body) => {
        if (pending.get(key) === promise) pending.delete(key);
        settle(body);
      };

      const originalJson = res.json.bind(res);
      res.json = (body) => {
        // error responses (e.g. from the error handler) must not be cached
        const ok = res.statusCode < 400;
        if (ok) {
          cache
            .set(key, JSON.stringify(body))
            .catch((error) =>
              req.app.locals.logger?.error(`Failed to cache ${key}: ${error}`)
            );
        }
        finish(ok ? body : undefined);
        return originalJson(body);
      };
      res.on('close', () => finish(undefined));
      next();
    } catch (err) {
      next(err);
    }
  };
}
