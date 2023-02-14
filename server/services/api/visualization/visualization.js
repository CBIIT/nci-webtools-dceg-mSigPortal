import { Router } from 'express';
import {
  getSeqmatrixData,
  getSeqmatrixOptions,
  getSeqmatrixSummary,
  getClusterData,
  getRefgenomeData,
} from '../../query.js';

async function querySeqmatrix(req, res, next) {
  try {
    const { limit, offset, rowMode, userId, ...query } = req.query;
    // const { study, cancer, strategy } = query;
    // if (!userId && (!study || !cancer || !strategy)) {
    //   throw Error(
    //     'Missing one or more of the following parameters: study, cancer, strategy'
    //   );
    // }
    const connection = userId
      ? req.app.locals.sqlite(userId, 'local')
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
      ? req.app.locals.sqlite(userId, 'local')
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
      ? req.app.locals.sqlite(userId, 'local')
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
    if (!userId) {
      throw 'This API is only available for user data calculations';
    }
    const connection = req.app.locals.sqlite(userId, 'local');
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
async function queryRefgenome(req, res, next) {
  try {
    const { limit, offset, rowMode, userId, ...query } = req.query;
    const connection = req.app.locals.connection;
    const columns = '*';
    const data = await getRefgenomeData(
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
router.get('/refgenome', queryRefgenome);

export { router, querySeqmatrix };
