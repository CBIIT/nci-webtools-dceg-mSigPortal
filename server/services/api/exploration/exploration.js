const { Router } = require('express');
const { randomUUID } = require('crypto');
const {
  getExposureData,
  getExposureOptions,
  getSignatureData,
  getSeqmatrixData,
} = require('../../query');
const { addBurden } = require('./burden');
const r = require('r-wrapper').async;
const path = require('path');
const fs = require('fs-extra');
const config = require('../../../config.json');
const logger = require('../../logger');

// config info for R functions
const rConfig = {
  s3Data: config.data.s3,
  bucket: config.data.bucket,
  localData: path.resolve(config.data.localData),
  wd: path.resolve(config.results.folder),
};

function alphaNumericSort(array) {
  return array.sort((a, b) => {
    return a.localeCompare(b, undefined, {
      numeric: true,
      sensitivity: 'base',
    });
  });
}
async function queryExposure(req, res, next) {
  try {
    const { limit, offset, orderByCluster, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = '*';
    const data = await getExposureData(
      connection,
      query,
      columns,
      limit,
      offset
    );

    if (orderByCluster) {
      try {
        const wrapper = await r('services/R/explorationWrapper.R', 'wrapper', {
          fn: 'hierarchicalClusterOrder',
          args: { data },
          config: { ...rConfig },
        });

        const { output: sampleClusterOrder } = JSON.parse(wrapper);

        const orderedData = data.sort(
          (a, b) =>
            sampleClusterOrder.indexOf(a.sample) -
            sampleClusterOrder.indexOf(b.sample)
        );

        res.json(addBurden(orderedData));
      } catch (err) {
        logger.error(`An error while calculating clusters`);
        next(err);
      }
    } else {
      res.json(addBurden(data));
    }
  } catch (error) {
    next(error);
  }
}

// query public exploration options for exploration tab
async function explorationOptions(req, res, next) {
  try {
    const { limit, offset, ...query } = req.query;
    const connection = req.app.locals.connection;

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

async function explorationSamples(req, res, next) {
  try {
    const { study, strategy, signatureSetName } = req.query;
    const connection = req.app.locals.connection;

    const query = { study, strategy, signatureSetName };
    const columns = ['sample', 'signatureName'];
    const data = await getExposureData(connection, query, columns);
    const samples = alphaNumericSort([...new Set(data.map((e) => e.sample))]);
    const signatureNames = alphaNumericSort([
      ...new Set(data.map((e) => e.signatureName)),
    ]);
    res.json({ samples, signatureNames });
  } catch (error) {
    next(error);
  }
}

async function explorationWrapper(req, res, next) {
  const { fn, args, projectID: id, type = 'calc' } = req.body;
  logger.debug('/explorationWrapper: ' + fn);
  // logger.debug('/explorationWrapper: %o', { ...req.body });

  const projectID = id ? id : type == 'calc' ? randomUUID() : false;
  // create directory for results if needed
  const savePath = projectID ? path.join(projectID, 'results', fn, '/') : null;
  if (projectID)
    fs.mkdirSync(path.join(rConfig.wd, savePath), { recursive: true });

  try {
    const wrapper = await r('services/R/explorationWrapper.R', 'wrapper', {
      fn,
      args,
      config: {
        ...rConfig,
        savePath,
        projectID,
      },
    });

    const { stdout, ...rest } = JSON.parse(wrapper);

    res.json({
      projectID,
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
  logger.debug('/msLandscape: ' + JSON.stringify(req.body));
  try {
    const {
      params_activity,
      params_spectrum,
      params_signature,
      projectID: id,
    } = req.body;

    const connection = req.app.locals.connection;
    const columns = '*';

    const exposureData = await getExposureData(
      connection,
      params_activity,
      columns
    );
    const signatureData = await getSignatureData(
      connection,
      params_signature,
      columns
    );
    const seqmatrixData = await getSeqmatrixData(
      connection,
      params_spectrum,
      columns
    );

    // console.log(exposureData[0]);
    // console.log(signatureData[0]);
    console.log(seqmatrixData[0]);
    const fn = 'msLandscape';
    const args = { exposureData, signatureData, seqmatrixData };
    const projectID = id || randomUUID();

    const savePath = path.join(projectID, 'results', fn, '/');
    fs.mkdirSync(path.join(rConfig.wd, savePath), { recursive: true });

    const wrapper = await r('services/R/explorationWrapper.R', 'wrapper', {
      fn,
      args,
      config: {
        ...rConfig,
        savePath,
        projectID,
      },
    });

    const { stdout, ...rest } = JSON.parse(wrapper);

    res.json({
      projectID,
      stdout,
      ...rest,
    });
  } catch (err) {
    logger.error(`/msLandscape: An error occured `);
    next(err);
  }
}

const router = Router();

router.get('/mutational_activity', queryExposure);
router.get('/mutational_activity_options', explorationOptions);
router.get('/explorationSamples', explorationSamples);
router.post('/explorationWrapper', explorationWrapper);
router.post('/msLandscape', msLandscape);

module.exports = { router, queryExposure };
