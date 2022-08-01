const { Router } = require('express');
const { getPublicationData } = require('../../query');
const { pickNonNullValues } = require('../../utils');

async function queryPublications(req, res, next) {
  try {
    const { limit, offset, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = '*';
    const data = await getPublicationData(
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

router.get('/publications', queryPublications);

module.exports = { router, queryPublications };
