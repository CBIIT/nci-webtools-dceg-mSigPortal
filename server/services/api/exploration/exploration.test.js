import request from 'supertest';
import express from 'express';
import knex from 'knex';
import { router } from './exploration.js';

const env = process.env;
const app = express();

app.use(express.json());
app.use(router);

app.locals.connection = knex({
  client: 'postgres',
  connection: {
    host: env.POSTGRES_HOST,
    port: env.POSTGRES_PORT,
    user: env.POSTGRES_USER,
    password: env.POSTGRES_PASS,
    database: env.POSTGRES_DB,
  },
});

describe('queryExposure', () => {
  test('should return an error if the required parameters userId, limit, offset, orderByCluster are not provided', async () => {
    const response = await request(app).get('/signature_activity').query({});

    expect(response.statusCode).toBe(400);
    expect(response.body).toEqual({
      error:
        'Missing one or more of the following parameters: userId, limit, offset, orderByCluster',
    });
  });

  test('should query the Exposure table with the provided parameters and return the results', async () => {
    const query = {
      study: 'PCAWG',
      strategy: 'WGS',
      signatureSetName: 'COSMIC_v3_Signatures_GRCh37_SBS96',
      limit: 10,
    };
    const response = await request(app).get('/signature_activity').query(query);

    expect(response.statusCode).toBe(200);
    expect(response.body.length).toBeGreaterThan(0);
  });
});

describe('exposureOptions', () => {
  it('should query the Exposure table with the provided parameters and return the results', async () => {
    const query = {
      study: 'PCAWG',
      limit: 10,
    };
    const response = await request(app)
      .get('/signature_activity_options')
      .query(query);

    expect(response.statusCode).toBe(200);
    expect(response.body.length).toBeGreaterThan(0);
  });
});
