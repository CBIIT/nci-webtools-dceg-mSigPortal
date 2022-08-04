const { Router } = require('express');
const { getPatternData } = require('../../query');

async function queryPattern(req, res, next) {
  try {
    const { limit, offset, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = '*';
    const data = await getPatternData(
      connection,
      query,
      columns,
      limit,
      offset
    );
    res.json(data);
  } catch (error) {
    next(error);
  }
}

const router = Router();

router.get('/pattern', queryPattern);

module.exports = { router, queryPattern };
