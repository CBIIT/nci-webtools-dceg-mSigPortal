const { Router } = require('express');
const { v4: uuidv4 } = require('uuid');
const { getSeqmatrixData, getSeqmatrixOptions } = require('../../query');

// query public seqmatrix data for visualization tab
async function visualizationOptions(req, res, next) {
  try {
    const { limit, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = '*';
    const data = await getSeqmatrixOptions(connection, query, columns, limit);
    const projectID = uuidv4();

    res.json(data);
  } catch (error) {
    next(error);
  }
}

async function visualizationSamples(req, res, next) {
  try {
    const { limit, ...query } = req.query;
    const connection = req.app.locals.connection;

    // const query = { study, strategy, cancer };
    const columns = ['profile', 'sample', 'profile', 'matrix'];
    const data = await getSeqmatrixData(connection, query, columns, limit);
    const projectID = uuidv4();

    res.json({ data, projectID });
  } catch (error) {
    next(error);
  }
}

async function querySeqmatrix(req, res, next) {
  try {
    const { limit, offset, rowMode, ...query } = req.query;
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

const router = Router();

router.get('/seqmatrix', querySeqmatrix);
router.get('/visualizationOptions', visualizationOptions);
router.get('/visualizationSamples', visualizationSamples);

module.exports = { router, querySeqmatrix };
