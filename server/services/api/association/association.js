const { Router } = require('express');
const { getAssociationData, getAssociationOptions } = require('../../query');

async function queryAssociation(req, res, next) {
  try {
    const { limit, offset, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = '*';
    const data = await getAssociationData(
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

async function associationOptions(req, res, next) {
  try {
    const { limit, offset, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = '*';
    const data = await getAssociationOptions(
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

const router = Router();

router.get('/association', queryAssociation);
router.get('/associationOptions', associationOptions);

module.exports = { router, queryAssociation, associationOptions };
