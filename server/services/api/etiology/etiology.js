const { Router } = require('express');
const { getEtiologyData } = require('../../query');
const { pickNonNullValues } = require('../../utils');

async function queryEtiology(req, res, next) {
  try {
    const { limit, offset, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = '*';
    const data = await getEtiologyData(
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

router.get('/etiology', queryEtiology);

module.exports = { router, queryEtiology };
