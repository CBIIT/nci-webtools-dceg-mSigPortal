const { Router } = require('express');

const {
  getSeqmatrixData,
  getSeqmatrixOptions,
  getSeqmatrixSummary,
  getClusterData,
} = require('../../query');

async function querySeqmatrix(req, res, next) {
  try {
    const { limit, offset, rowMode, userId, ...query } = req.query;
    const { study, cancer, strategy } = query;
    if (!userId && (!study || !cancer || !strategy)) {
      throw 'Missing one or more of the following parameters: study, cancer, strategy';
    }

    const connection = userId
      ? req.app.locals.sqlite(userId, 'visualization')
      : req.app.locals.connection;

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
    const { limit, offset, userId, ...query } = req.query;

    const connection = userId
      ? req.app.locals.sqlite(userId, 'visualization')
      : req.app.locals.connection;

    const columns = '*';
    const data = await getSeqmatrixOptions(
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

async function seqmatrixSummary(req, res, next) {
  try {
    const { limit, offset, rowMode, userId, ...query } = req.query;
    const { study, cancer, strategy } = query;
    if (!userId && (!study || !cancer || !strategy)) {
      throw 'Missing one or more of the following parameters: study, cancer, strategy';
    }

    const connection = userId
      ? req.app.locals.sqlite(userId, 'visualization')
      : req.app.locals.connection;

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

async function queryCluster(req, res, next) {
  try {
    const { limit, offset, rowMode, userId, ...query } = req.query;
    const { study, cancer, strategy } = query;
    if (!userId && (!study || !cancer || !strategy)) {
      throw 'Missing one or more of the following parameters: study, cancer, strategy';
    }

    const connection = req.app.locals.sqlite(userId, 'visualization');

    const columns = '*';
    const data = await getClusterData(
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

router.get('/mutational_spectrum', querySeqmatrix);
router.get('/mutational_spectrum_options', seqmatrixOptions);
router.get('/mutational_spectrum_summary', seqmatrixSummary);
router.get('/cluster', queryCluster);

module.exports = { router, querySeqmatrix };
