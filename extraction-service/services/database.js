import { createReadStream, existsSync } from "fs";
import { format } from "util";
import { parse } from "csv-parse";
import knex from "knex";

export async function createDatabaseFromFiles(tables, databaseFile) {
  const connection = getSqliteConnection(databaseFile);
  for (const { name, file } of tables) {
    if (existsSync(file)) {
      const delimiter = ["\t", ...createDelimiters(" ", 100)];
      await createSqliteTableFromFile(connection, name, file, { delimiter });
    }
  }
}

export async function createSqliteTableFromFile(connection, table, filepath, parseOptions = {}) {
  const parser = getDelimitedFileParser(filepath, parseOptions);
  const header = await getHeader(filepath, parseOptions);
  const tableDef = getTableDef(header);
  return await createSqliteTable(connection, table, tableDef, parser);
}

export function getSqliteConnection(filepath, options = {}) {
  return knex({
    client: "better-sqlite3",
    useNullAsDefault: true,
    connection: {
      filename: filepath,
      ...options,
    },
  });
}

export function createDelimiters(character, maxLength) {
  const delimiters = [];
  for (let i = maxLength; i > 0; i--) {
    const delimiter = character.repeat(i);
    delimiters.push(delimiter);
  }
  return delimiters;
}

export function getDelimitedFileParser(filepath, options = {}) {
  const parseOptions = {
    delimiter: createDelimiters(" ", 100),
    trim: true,
    cast: true,
    cast_date: true,
    columns: true,
    skip_empty_lines: true,
    ...options,
  };
  const readStream = createReadStream(filepath);
  return readStream.pipe(parse(parseOptions));
}

export async function createSqliteTable(connection, table, tableDef, iterator, bufferSize = 500) {
  await connection.schema.dropTableIfExists(table);
  await connection.schema.createTable(table, tableDef);
  return await importTable(connection, table, iterator, bufferSize);
}

async function peek(iterator) {
  const next = iterator.next();
  const resetIterator = function* () {
    if (next.done) return;
    yield next.value;
    yield* iterator;
  };
  return { value: await next.value, resetIterator };
}

export async function importTable(connection, table, iterator, bufferSize = 500, logger = console) {
  let count = 0;
  let buffer = [];

  async function flushBuffer() {
    try {
      await connection.batchInsert(table, buffer);
      count += buffer.length;
      buffer = [];
      logger.info(`Imported ${count} rows`);
    } catch (error) {
      // batchInsert exceptions do not return the specific records that failed
      // so we need to check each record individually
      for (let record of buffer) {
        try {
          await connection(table).insert(record);
        } catch (error) {
          logger.error(record);
          throw error;
        }
      }
      throw error;
    }
  }

  for await (const record of iterator) {
    buffer.push(record);
    if (buffer.length >= bufferSize) {
      await flushBuffer();
    }
  }

  await flushBuffer();
  return count;
}

export async function getHeader(filepath, parseOptions = {}) {
  const parser = getDelimitedFileParser(filepath, { ...parseOptions, to_line: 2 });
  const header = await parser.iterator().next();
  return await header.value;
}

export function getColumnDefs(header) {
  return Object.keys(header).map((key) => ({
    name: key,
    type: getType(header[key]),
  }));
}

export function getTableDef(header) {
  const columnDefs = getColumnDefs(header);
  return (table) => columnDefs.forEach((c) => table[c.type](c.name));
}

export function getType(value) {
  if (typeof value === "number") {
    return "double";
  } else if (typeof value === "string" || value === undefined) {
    return "text";
  } else if (typeof value === "boolean") {
    return "boolean";
  } else if (value instanceof Date) {
    return "datetime";
  }
  throw new Error(`Unsupported type for value: ${format(value)}`);
}
