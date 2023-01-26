import { pickBy } from 'lodash-es';

export function getData(
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

export function getSignatureData(
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

export function getQueryMiddleware(queryFunction) {
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
