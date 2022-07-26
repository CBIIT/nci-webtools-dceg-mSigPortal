const { Router } = require('express');
const { getAssociationData } = require('../../query');

async function queryAssociation(req, res, next) {
  try {
    const { limit, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = '*';
    const data = await getAssociationData(connection, query, columns, limit);
    res.json(data);
  } catch (error) {
    next(error);
  }
}

const router = Router();

router.get('/association', queryAssociation);

module.exports = { router, queryAssociation };
