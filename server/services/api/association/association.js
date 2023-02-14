import { Router } from 'express';
import { getAssociationData, getAssociationOptions } from '../../query.js';

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
router.get('/signature_association', queryAssociation);
router.get('/signature_association_options', associationOptions);

export { router, queryAssociation, associationOptions };
