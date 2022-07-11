const { pickBy } = require('lodash');

function getData(
  connection,
  table,
  query,
  columns = '*',
  limit = 1000000,
  offset = 0
) {
  const conditions = pickBy(query, (v) => v !== undefined);
  return connection
    .distinct(columns)
    .from(table)
    .where(conditions)
    .limit(limit)
    .offset(offset);
}

function getAssociationData(
  connection,
  query,
  columns = '*',
  limit = 100000,
  offset = 0
) {
  return getData(connection, 'association', query, columns, limit, offset);
}

function getExposureData(
  connection,
  query,
  columns = '*',
  limit = 100000,
  offset = 0
) {
  return getData(connection, 'exposure', query, columns, limit, offset);
}

function getSeqmatrixData(
  connection,
  query,
  columns = '*',
  limit = 100000,
  offset = 0
) {
  return getData(connection, 'seqmatrix', query, columns, limit, offset);
}

function getSignatureData(
  connection,
  query,
  columns = '*',
  limit = 100000,
  offset = 0
) {
  return getData(connection, 'signature', query, columns, limit, offset);
}

module.exports = {
  getData,
  getAssociationData,
  getExposureData,
  getSeqmatrixData,
  getSignatureData,
};
