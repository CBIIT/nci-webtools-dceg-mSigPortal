const { Router } = require('express');
const { getSignatureData, getSignatureOptions } = require('../../query');

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

// query public exploration options for exploration tab
async function signatureOptions(req, res, next) {
  try {
    const { study, strategy, signatureSetName, cancer } = req.query;
    const connection = req.app.locals.connection;

    const query = { study, strategy, signatureSetName, cancer };
    const columns = ['study', 'strategy', 'cancer', 'signatureSetName'];
    const data = await getSignatureOptions(connection, query, columns);
    res.json(data);
  } catch (error) {
    next(error);
  }
}

const router = Router();

router.get('/signature', querySignature);
router.get('/signatureOptions', signatureOptions);

module.exports = { router, querySignature };
