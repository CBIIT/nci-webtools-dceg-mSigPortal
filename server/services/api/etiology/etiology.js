import { Router } from 'express';
import {
  getEtiologyOptions,
  getEtiologyData,
  getEtiologyOrganData,
} from '../../query.js';
import { pickNonNullValues } from '../../utils.js';

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

async function queryEtiologyOrgan(req, res, next) {
  try {
    const { limit, offset, ...query } = req.query;
    const connection = req.app.locals.connection;
    const columns = '*';
    const data = await getEtiologyOrganData(
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
router.get('/signature_etiology_options', queryEtiologyOptions);
router.get('/signature_etiology', queryEtiology);
router.get('/signature_etiology_organ', queryEtiologyOrgan);

export { router, queryEtiologyOptions };
