const { Router } = require('express');
const {
  getEtiologyOptions,
  getEtiologyData,
  getEtiologyOrganData,
  getEtiologySignatureData,
} = require('../../query');
const { pickNonNullValues } = require('../../utils');

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

async function queryEtiologySignature(req, res, next) {
  try {
    const { limit, offset, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = '*';
    const data = await getEtiologySignatureData(
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

router.get('/signature_etiology_signature', queryEtiologySignature);

module.exports = { router, queryEtiologyOptions };
