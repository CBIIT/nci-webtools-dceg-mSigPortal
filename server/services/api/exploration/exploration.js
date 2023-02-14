import express from 'express';
import { randomUUID } from 'crypto';
import {
  getExposureData,
  getExposureOptions,
  getSignatureData,
  getSeqmatrixData,
} from '../../query.js';
import { addBurden } from './burden.js';
import rWrapper from 'r-wrapper';
import path from 'path';
import fs from 'fs-extra';
import config from '../../../config.json' assert { type: 'json' };
import logger from '../../logger.js';
const { Router } = express;
const r = rWrapper.async;

// config info for R functions
const rConfig = {
  s3Data: config.data.s3,
  bucket: config.data.bucket,
  localData: path.resolve(config.data.localData),
  wd: path.resolve(config.results.folder),
};
async function queryExposure(req, res, next) {
  try {
    const { userId, limit, offset, orderByCluster, ...query } = req.query;
    const connection = userId
      ? req.app.locals.sqlite(userId, 'local')
      : req.app.locals.connection;
    const columns = '*';
    const data = await getExposureData(
      connection,
      query,
      columns,
      limit,
      offset
    );
    res.json(addBurden(data));
  } catch (error) {
    next(error);
  }
}
// query exposure options for exploration tab
async function exposureOptions(req, res, next) {
  try {
    const { userId, limit, offset, ...query } = req.query;
    const connection = userId
      ? req.app.locals.sqlite(userId, 'local')
      : req.app.locals.connection;
    const columns = ['study', 'strategy', 'cancer', 'signatureSetName'];
    const data = await getExposureOptions(
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
async function explorationWrapper(req, res, next) {
  const { fn, args, id, type = 'calc' } = req.body;
  logger.debug('/explorationWrapper: ' + fn);
  // logger.debug('/explorationWrapper: %o', { ...req.body });
  const sessionId = id ? id : type == 'calc' ? randomUUID() : false;
  // create directory for results if needed
  const savePath = id ? path.join(id, 'results', fn, '/') : null;
  if (sessionId)
    fs.mkdirSync(path.join(rConfig.wd, savePath), { recursive: true });
  try {
    const wrapper = await r('services/R/explorationWrapper.R', 'wrapper', {
      fn,
      args,
      config: {
        ...rConfig,
        savePath,
        id,
      },
    });
    const { stdout, ...rest } = JSON.parse(wrapper);
    res.json({
      id,
      stdout,
      ...rest,
    });
  } catch (err) {
    logger.error(`/explorationCalc: An error occured with fn: ${fn}`);
    next(err);
  }
}
// query 3 separate tables and return exposure data sorted by clusters and cosine similarity
async function msLandscape(req, res, next) {
  try {
    const { study, strategy, signatureSetName, cancer, userId } = req.query;
    const connection = userId
      ? req.app.locals.sqlite(userId, 'local')
      : req.app.locals.connection;
    const columns = '*';
    const limit = false;
    const exposureData = await getExposureData(
      connection,
      { study, strategy, signatureSetName, cancer },
      columns,
      limit
    );
    const signatureData = await getSignatureData(
      connection,
      { signatureSetName },
      columns,
      limit
    );
    const seqmatrixData = await getSeqmatrixData(
      connection,
      { study, strategy, cancer },
      columns,
      limit
    );
    const fn = 'msLandscape';
    const args = { exposureData, signatureData, seqmatrixData };
    const id = userId || randomUUID();
    const wrapper = await r('services/R/explorationWrapper.R', 'wrapper', {
      fn,
      args,
      config: {
        ...rConfig,
        id,
      },
    });
    const { stdout, ...rest } = JSON.parse(wrapper);
    res.json({
      id,
      stdout,
      ...rest,
    });
  } catch (err) {
    logger.error(`/msLandscape: An error occured `);
    next(err);
  }
}
async function msDecomposition(req, res, next) {
  try {
    const { study, strategy, signatureSetName, cancer, userId } = req.query;
    const connection = userId
      ? req.app.locals.sqlite(userId, 'local')
      : req.app.locals.connection;
    const columns = '*';
    const limit = false;
    const exposureData = await getExposureData(
      connection,
      { study, strategy, signatureSetName, cancer },
      columns,
      limit
    );
    const signatureData = await getSignatureData(
      connection,
      { signatureSetName },
      columns,
      limit
    );
    const seqmatrixData = await getSeqmatrixData(
      connection,
      { study, strategy, cancer },
      columns,
      limit
    );
    const fn = 'msDecomposition';
    const args = { exposureData, signatureData, seqmatrixData };
    const id = userId || randomUUID();
    const wrapper = await r('services/R/explorationWrapper.R', 'wrapper', {
      fn,
      args,
      config: {
        ...rConfig,
        id,
      },
    });
    const { stdout, ...rest } = JSON.parse(wrapper);
    res.json({
      id,
      stdout,
      ...rest,
    });
  } catch (err) {
    logger.error(`/msDecomposition: An error occured `);
    next(err);
  }
}
// query signature data and calculate cosine similarity
async function cosineSimilarity(req, res, next) {
  try {
    const { signatureSetName, userId, ...params } = req.query;
    const connection = userId
      ? req.app.locals.sqlite(userId, 'local')
      : req.app.locals.connection;
    const columns = '*';
    const limit = false;
    const signatureData1 = await getSignatureData(
      connection,
      { ...params, signatureSetName: signatureSetName.split(';')[0] },
      columns,
      limit
    );
    const signatureData2 = await getSignatureData(
      connection,
      { ...params, signatureSetName: signatureSetName.split(';')[1] },
      columns,
      limit
    );
    const fn = 'cosineSimilarity';
    const args = { signatureData1, signatureData2 };
    const wrapper = await r('services/R/explorationWrapper.R', 'wrapper', {
      fn,
      args,
      config: {
        ...rConfig,
      },
    });
    const { stdout, ...rest } = JSON.parse(wrapper);
    res.json({ stdout, ...rest });
  } catch (err) {
    logger.error(`/msLandscape: An error occured `);
    next(err);
  }
}
const router = Router();
router.get('/signature_activity', queryExposure);
router.get('/signature_activity_options', exposureOptions);
router.post('/explorationWrapper', explorationWrapper);
router.get('/signature_landscape', msLandscape);
router.get('/signature_decomposition', msDecomposition);
router.get('/signature_cosine_similarity', cosineSimilarity);
export { router };
export { queryExposure };
export default {
  router,
  queryExposure,
};
