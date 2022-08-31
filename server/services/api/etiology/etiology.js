const { Router } = require('express');
const { getEtiologyOptions } = require('../../query');
const { pickNonNullValues } = require('../../utils');

async function queryEtiologyOptions(req, res, next) {
  try {
    const { limit, offset, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = '*';
    const data = await getEtiologyOptions(
      connection,
      query,
      columns,
      limit,
      offset
    );
    const records = data.map(pickNonNullValues);
    res.json(records);
  } catch (error) {
    next(error);
  }
}

const router = Router();

router.get('/etiologyOptions', queryEtiologyOptions);

module.exports = { router, queryEtiologyOptions };
