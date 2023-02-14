import { pickBy, cond } from 'lodash-es';

function getData(
  connection,
  table,
  query,
  columns = '*',
  limit,
  offset = 0,
  rowMode = 'object',
  distinct = false
) {
  const conditions = pickBy(
    query,
    (v) => v && !v.includes('%') && !v.includes('*ALL')
  );
  const patterns = pickBy(query, (v) => v && v.includes('%'));

  let sqlQuery = connection
    .select(columns)
    .from(table)
    // .where(conditions)
    .offset(offset, rowMode)
    .options({ rowMode: rowMode });

  if (limit) {
    sqlQuery = sqlQuery.limit(limit || 100000);
  }
  if (distinct) {
    sqlQuery = sqlQuery.distinct(columns);
  }

  // apply where conditions to query
  // use WHERE IN query on conditions delimited by semi-colons (;)
  if (conditions) {
    Object.entries(conditions).forEach(([column, values]) => {
      const splitValues = values.split(';');
      if (splitValues.length > 1) {
        sqlQuery.whereIn(
          column,
          splitValues.map((e) => e.trim())
        );
      } else {
        splitValues.forEach((v) => {
          sqlQuery = sqlQuery.andWhere(column, v.trim());
        });
      }
    });
  }

  if (patterns) {
    Object.entries(patterns).forEach(([column, values]) => {
      sqlQuery = sqlQuery.andWhere(column, 'like', values.trim());
    });
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

function getAssociationOptions(
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
    'associationOption',
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

function getSignatureSummary(
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
    'signatureSummary',
    query,
    columns,
    limit,
    offset,
    rowMode,
    distinct
  );
}

function getEtiologyOptions(
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
    'etiologyOptions',
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

function getEtiologyOrganData(
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
    'etiologyOrgan',
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

function getPatternData(
  connection,
  query,
  columns = '*',
  limit = 400000,
  offset = 0,
  rowMode = 'object',
  distinct = false
) {
  return getData(
    connection,
    'pattern',
    query,
    columns,
    limit,
    offset,
    rowMode,
    distinct
  );
}

function getClusterData(
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
    'cluster',
    query,
    columns,
    limit,
    offset,
    rowMode,
    distinct
  );
}

function getRefgenomeData(
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
    'refgenome',
    query,
    columns,
    limit,
    offset,
    rowMode,
    distinct
  );
}

function getQueryMiddleware(queryFunction) {
  return async function (req, res, next) {
    try {
      const connection = req.app.locals.connection;
      let {
        query,
        columns = '*',
        limit = 200000,
        offset = 0,
        rowMode = 'object',
        distinct = 'false',
      } = req.body;
      const response = await queryFunction(
        connection,
        query,
        columns,
        limit,
        offset,
        rowMode,
        distinct
      );
      console.log(response);
      res.json(response);
    } catch (e) {
      next(e);
    }
  };
}

export {
  getData,
  getAssociationData,
  getAssociationOptions,
  getExposureData,
  getExposureOptions,
  getSeqmatrixData,
  getSeqmatrixOptions,
  getSeqmatrixSummary,
  getSignatureData,
  getSignatureOptions,
  getSignatureSummary,
  getEtiologyOptions,
  getEtiologyData,
  getEtiologyOrganData,
  getPublicationData,
  getPatternData,
  getClusterData,
  getRefgenomeData,
  getQueryMiddleware,
};
