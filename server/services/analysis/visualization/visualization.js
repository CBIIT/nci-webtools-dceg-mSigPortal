const { v4: uuidv4 } = require('uuid');
const { getSeqmatrixData, getSignatureData } = require('../../query');
const { SBS, DBS, ID } = require('./mutationalProfiles');

// query public seqmatrix data for visualization tab
async function visualizationOptions(req, res, next) {
  try {
    const { limit, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = ['study', 'cancer', 'strategy'];
    const data = await getSeqmatrixData(connection, query, columns, limit);
    const projectID = uuidv4();

    res.json({ data, projectID });
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

async function mutationalProfiles(req, res, next) {
  try {
    const { limit, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = ['mutationType', 'mutations'];
    const data = await getSeqmatrixData(connection, query, columns, limit);

    if (query.profile == 'SBS') {
      res.json(SBS(data));
    } else if (query.profile == 'DBS') {
      res.json(data);
    } else if (query.profile == 'ID') {
      res.json(data);
    } else {
      throw 'mutationalProfiles: unsupported profile';
    }
  } catch (error) {
    next(error);
  }
}

async function profilerSummary(req, res, next) {
  try {
    const { limit, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = ['sample', 'profile', 'mutations'];
    const data = await getSeqmatrixData(connection, query, columns, limit);
    res.json(data);
  } catch (error) {
    next(error);
  }
}

async function querySeqmatrix(req, res, next) {
  try {
    const { limit, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = ['mutationType', 'mutations'];
    const data = await getSeqmatrixData(connection, query, columns, limit);
    res.json(data);
  } catch (error) {
    next(error);
  }
}

async function querySignature(req, res, next) {
  try {
    const { limit, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = !query.signatureSetName
      ? ['signatureSetName']
      : [
          'strandInfo',
          'strand',
          'signatureName',
          'mutationType',
          'contribution',
        ];
    const data = await getSignatureData(connection, query, columns, limit);
    res.json(data);
  } catch (error) {
    next(error);
  }
}

module.exports = {
  visualizationOptions,
  visualizationSamples,
  mutationalProfiles,
  profilerSummary,
  querySeqmatrix,
  querySignature,
};