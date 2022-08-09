const { Router } = require('express');
const { v4: uuidv4 } = require('uuid');
const {
  getSeqmatrixData,
  getSeqmatrixOptions,
  getSeqmatrixSummary,
} = require('../../query');

async function querySeqmatrix(req, res, next) {
  try {
    const { limit, offset, rowMode, ...query } = req.query;
    const { study, cancer, strategy } = query;
    if (!study || !cancer || !strategy) {
      throw 'Missing one or more of the following parameters: study, cancer, strategy';
    }

    const connection = req.app.locals.connection;
    const columns = '*';
    const data = await getSeqmatrixData(
      connection,
      query,
      columns,
      limit,
      offset,
      rowMode
    );
    res.json(data);
  } catch (error) {
    next(error);
  }
}

// query public seqmatrix data for visualization tab
async function seqmatrixOptions(req, res, next) {
  try {
    const { limit, offset, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = '*';
    const data = await getSeqmatrixOptions(connection, query, columns, limit, offset);
    const projectID = uuidv4();

    res.json(data);
  } catch (error) {
    next(error);
  }
}

async function seqmatrixSummary(req, res, next) {
  try {
    const { limit, offset, rowMode, ...query } = req.query;
    const { study, cancer, strategy } = query;
    if (!study || !cancer || !strategy) {
      throw 'Missing one or more of the following parameters: study, cancer, strategy';
    }

    const connection = req.app.locals.connection;
    const columns = '*';
    const data = await getSeqmatrixSummary(
      connection,
      query,
      columns,
      limit,
      offset,
      rowMode
    );
    res.json(data);
  } catch (error) {
    next(error);
  }
}

const router = Router();

router.get('/seqmatrix', querySeqmatrix);
router.get('/seqmatrixOptions', seqmatrixOptions);
router.get('/seqmatrixSummary', seqmatrixSummary);

module.exports = { router, querySeqmatrix };
