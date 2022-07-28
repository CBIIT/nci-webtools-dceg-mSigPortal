const { pickBy } = require('lodash');

function getData(
  connection,
  table,
  query,
  columns = '*',
  limit = 1000000,
  offset = 0,
  rowMode = 'object'
) {
  const conditions = pickBy(query, (v) => v !== undefined);
  return connection
    .distinct(columns)
    .from(table)
    .where(conditions)
    .limit(limit)
    .offset(offset, rowMode)
    .options({ rowMode: rowMode });
}

function getAssociationData(
  connection,
  query,
  columns = '*',
  limit = 200000,
  offset = 0,
  rowMode = 'object'
) {
  return getData(
    connection,
    'association',
    query,
    columns,
    limit,
    offset,
    rowMode
  );
}

function getExposureData(
  connection,
  query,
  columns = '*',
  limit = 200000,
  offset = 0,
  rowMode = 'object'
) {
  return getData(
    connection,
    'exposure',
    query,
    columns,
    limit,
    offset,
    rowMode
  );
}

function getSeqmatrixData(
  connection,
  query,
  columns = '*',
  limit = 200000,
  offset = 0,
  rowMode = 'object'
) {
  return getData(
    connection,
    'seqmatrix',
    query,
    columns,
    limit,
    offset,
    rowMode
  );
}

function getSeqmatrixOptions(
  connection,
  query,
  columns = '*',
  limit = 200000,
  offset = 0,
  rowMode = 'object'
) {
  return getData(
    connection,
    'seqmatrixOption',
    query,
    columns,
    limit,
    offset,
    rowMode
  );
}

function getSignatureData(
  connection,
  query,
  columns = '*',
  limit = 200000,
  offset = 0,
  rowMode = 'object'
) {
  return getData(
    connection,
    'signature',
    query,
    columns,
    limit,
    offset,
    rowMode
  );
}

module.exports = {
  getData,
  getAssociationData,
  getExposureData,
  getSeqmatrixData,
  getSeqmatrixOptions,
  getSignatureData,
};
