import knex from "knex";
import postgres from "pg";
import copyStreams from "pg-copy-streams";
import { S3Client } from "@aws-sdk/client-s3";
import { LocalProvider } from "./providers/localProvider.js";
import { S3Provider } from "./providers/s3Provider.js";

export function createConnection(
  args = {
    host: "localhost",
    port: 5432,
    user: "msigportal",
    password: "msigportal",
    database: "msigportal",
  }
) {
  return knex({
    client: "pg",
    connection: args,
    pool: {
      min: 2,
      max: 20,
      reapIntervalMillis: 50,
      acquireTimeoutMillis: 100 * 1000,
      createTimeoutMillis: 100 * 1000,
      propagateCreateError: false,
    },
  });
}

export async function createPostgresConnection(
  args = {
    host: "localhost",
    port: 5432,
    user: "msigportal",
    password: "msigportal",
    database: "msigportal",
  }
) {
  const client = new postgres.Client(args);
  await client.connect();
  return client;
}

export function loadAwsCredentials(config) {
  const { region, accessKeyId, secretAccessKey } = config;

  if (region) {
    process.env.AWS_REGION = region;
    process.env.AWS_DEFAULT_REGION = region;
  }

  if (accessKeyId) {
    process.env.AWS_ACCESS_KEY_ID = accessKeyId;
  }

  if (secretAccessKey) {
    process.env.AWS_SECRET_ACCESS_KEY = secretAccessKey;
  }
}

export async function recreateTable(connection, name, schemaFn) {
  await connection.schema.dropTableIfExists(name);
  await connection.schema.createTable(name, (table) => schemaFn(table, connection));
  return true;
}

export async function initializeSchema(connectionConfig, schema) {
  const connection = createConnection(connectionConfig);
  const tables = schema.filter(({ type }) => !type || type === "table");
  const materializedViews = schema.filter(({ type }) => type === "materializedView");
  const indexedTables = schema.filter(s => typeof s.index === 'function')

  // drop tables in reverse order to avoid foreign key constraints
  for (const { name } of [...materializedViews.reverse()]) {
    await connection.schema.dropMaterializedViewIfExists(name);
  }

  for (const { name } of [...tables].reverse()) {
    await connection.schema.dropTableIfExists(name);
  }

  for (const { name, schema } of tables) {
    await connection.schema.createTable(name, (table) => schema(table, connection));
  }

  for (const { name, schema } of materializedViews) {
    await connection.schema.createMaterializedView(name, (view) => schema(view, connection));
  }

  for (const { name, index } of indexedTables) {
    await connection.schema.table(name, index);
  }

  return true;
}

export async function initializeSchemaForImport(connection, schema, sources) {
  const shouldRecreateTable = (table) => sources.find((s) => s.table === table.name || table.dependsOn?.includes(s.table));
  const importSchema = schema.filter(shouldRecreateTable);
  return await initializeSchema(connection, importSchema);
}

export function importPostgresTable(connection, inputStream, table, columns) {
  return new Promise((resolve, reject) => {
    const columnSql = columns.map(c => `"${c}"`).join(',');
    const copyStream = copyStreams.from(`COPY "${table}" (${columnSql}) FROM STDIN DELIMITER ',' CSV HEADER`)
    const stream = connection.query(copyStream)
    inputStream.on('error', reject);
    stream.on('error', reject);
    stream.on('finish', () => resolve(copyStream.rowCount));
    inputStream.pipe(stream)
  });
}

export async function withDuration(fn) {
  const start = Date.now();
  const results = await fn();
  const end = Date.now();
  const duration = (end - start) / 1000;
  return { results, duration };
}

export function getSourceProvider(providerName, providerArgs) {
  switch (providerName) {
    case "local":
      return new LocalProvider(...providerArgs);
    case "s3":
      return new S3Provider(new S3Client(), ...providerArgs);
    default:
      throw new Error(`Unknown provider: ${providerName}`);
  }
}