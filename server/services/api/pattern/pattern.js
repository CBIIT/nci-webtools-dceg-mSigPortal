import { Router } from 'express';
import logger from '../../logger.js';
import { getPatternData } from '../../query.js';
import { wrapper } from '../visualization/userVisualization.js';
import path from 'path';
import config from '../../../config.json' assert { type: 'json' };

// for public data, query the pattern table
// for user data, generate pattern data in R
async function queryPattern(req, res, next) {
  try {
    const { limit, offset, proportion, ...query } = req.query;
    if (query.study) {
      // public data query
      const connection = req.app.locals.connection;
      const columns = '*';
      let sql = getPatternData(connection, query, columns, limit, offset);
      if (proportion) sql = sql.where('n1', '>', proportion);
      res.json(await sql);
    } else {
      // generate for user data
      const params = {
        fn: 'mpeaUser',
        args: {
          pattern: query.pattern,
          proportion: proportion,
          matrixFile: path.join(
            config.results.folder,
            query.userId,
            query.matrixFile
          ),
        },
        config: {},
      };
      const { output, ...logs } = await wrapper('wrapper', params);
      logger.debug(logs);
      res.json(output);
    }
  } catch (error) {
    next(error);
  }
}

const router = Router();
router.get('/pattern', queryPattern);

export { router, queryPattern };
