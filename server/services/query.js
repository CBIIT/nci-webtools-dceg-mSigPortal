const { pickBy } = require('lodash');

function getData(
  connection,
  table,
  query,
  columns = '*',
  limit = 1000000,
  offset = 0,
  rowMode = 'object',
  distinct = false
) {
  const conditions = pickBy(query, (v) => v !== undefined);
  let sqlQuery = connection
    .select(columns)
    .from(table)
    .where(conditions)
    .limit(limit)
    .offset(offset, rowMode)
    .options({ rowMode: rowMode });

  if (distinct) {
    sqlQuery = sqlQuery.distinct(columns);
  }

  return sqlQuery;
}

function getAssociationData(
  connection,
  query,
  columns = '*',
  limit = 200000,
  offset = 0,
  rowMode = 'object',
  distinct = false
) {
  return getData(
    connection,
    'association',
    query,
    columns,
    limit,
    offset,
    rowMode,
    distinct
  );
}

function getExposureData(
  connection,
  query,
  columns = '*',
  limit = 200000,
  offset = 0,
  rowMode = 'object',
  distinct = false
) {
  return getData(
    connection,
    'exposure',
    query,
    columns,
    limit,
    offset,
    rowMode,
    distinct
  );
}

function getExposureOptions(
  connection,
  query,
  columns = '*',
  limit = 200000,
  offset = 0,
  rowMode = 'object',
  distinct = false
) {
  return getData(
    connection,
    'exposureOption',
    query,
    columns,
    limit,
    offset,
    rowMode,
    distinct
  );
}

function getSeqmatrixData(
  connection,
  query,
  columns = '*',
  limit = 200000,
  offset = 0,
  rowMode = 'object',
  distinct = false
) {
  return getData(
    connection,
    'seqmatrix',
    query,
    columns,
    limit,
    offset,
    rowMode,
    distinct
  );
}

function getSeqmatrixOptions(
  connection,
  query,
  columns = '*',
  limit = 200000,
  offset = 0,
  rowMode = 'object',
  distinct = false
) {
  return getData(
    connection,
    'seqmatrixOption',
    query,
    columns,
    limit,
    offset,
    rowMode,
    distinct
  );
}

function getSeqmatrixSummary(
  connection,
  query,
  columns = '*',
  limit = 200000,
  offset = 0,
  rowMode = 'object',
  distinct = false
) {
  return getData(
    connection,
    'seqmatrixSummary',
    query,
    columns,
    limit,
    offset,
    rowMode,
    distinct
  );
}

function getSignatureData(
  connection,
  query,
  columns = '*',
  limit = 200000,
  offset = 0,
  rowMode = 'object',
  distinct = false
) {
  return getData(
    connection,
    'signature',
    query,
    columns,
    limit,
    offset,
    rowMode,
    distinct
  );
}

function getSignatureOptions(
  connection,
  query,
  columns = '*',
  limit = 200000,
  offset = 0,
  rowMode = 'object',
  distinct = false
) {
  return getData(
    connection,
    'signatureOption',
    query,
    columns,
    limit,
    offset,
    rowMode,
    distinct
  );
}

function getEtiologyData(
  connection,
  query,
  columns = '*',
  limit = 200000,
  offset = 0,
  rowMode = 'object',
  distinct = false
) {
  return getData(
    connection,
    'etiology',
    query,
    columns,
    limit,
    offset,
    rowMode,
    distinct
  );
}

function getPublicationData(
  connection,
  query,
  columns = '*',
  limit = 200000,
  offset = 0,
  rowMode = 'object',
  distinct = false
) {
  return getData(
    connection,
    'publication',
    query,
    columns,
    limit,
    offset,
    rowMode,
    distinct
  );
}

module.exports = {
  getData,
  getAssociationData,
  getExposureData,
  getExposureOptions,
  getSeqmatrixData,
  getSeqmatrixOptions,
  getSeqmatrixSummary,
  getSignatureData,
  getSignatureOptions,
  getEtiologyData,
  getPublicationData,
};
