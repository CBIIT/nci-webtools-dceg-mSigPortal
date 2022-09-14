const { Router } = require('express');
const {
  getEtiologyOptions,
  getEtiologyData,
  getEtiologyOrganData,
} = require('../../query');
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

async function queryEtiology(req, res, next) {
  try {
    const { limit, offset, type, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = '*';
    const data =
      type == 'organ'
        ? await getEtiologyOrganData(connection, query, columns, limit, offset)
        : await getEtiologyData(connection, query, columns, limit, offset);
    const records = data.map(pickNonNullValues);
    res.json(records);
  } catch (error) {
    next(error);
  }
}

const router = Router();

router.get('/signatureEtiologyOptions', queryEtiologyOptions);

router.get('/signatureEtiology', queryEtiology);

module.exports = { router, queryEtiologyOptions };
