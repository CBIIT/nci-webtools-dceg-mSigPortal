const { Router } = require('express');
const { getPatternData } = require('../../query');

async function queryPattern(req, res, next) {
  try {
    const { limit, offset, proportion, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = '*';
    let sql = getPatternData(connection, query, columns, limit, offset);
    if (proportion) sql = sql.where('n1', '>', proportion);

    res.json(await sql);
  } catch (error) {
    next(error);
  }
}

const router = Router();

router.get('/pattern', queryPattern);

module.exports = { router, queryPattern };
