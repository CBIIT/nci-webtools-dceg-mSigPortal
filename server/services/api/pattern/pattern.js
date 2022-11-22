const { Router } = require('express');
const logger = require('../../logger');
const { getPatternData } = require('../../query');
const { wrapper } = require('../visualization/userVisualization');
const path = require('path');
const config = require('../../../config.json');

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

module.exports = { router, queryPattern };
