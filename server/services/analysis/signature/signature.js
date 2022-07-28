const { Router } = require('express');
const { v4: uuidv4 } = require('uuid');
const { getSignatureData } = require('../../query');

async function querySignature(req, res, next) {
  try {
    const { limit, offset, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = '*';
    const rowMode = 'object';
    const distinct = true;
    const data = await getSignatureData(
      connection,
      query,
      columns,
      limit,
      offset,
      rowMode,
      distinct
    );
    res.json(data);
  } catch (error) {
    next(error);
  }
}

const router = Router();

router.get('/signature', querySignature);

module.exports = { router, querySignature };
