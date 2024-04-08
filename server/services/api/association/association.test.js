import request from 'supertest';
import express from 'express';
import knex from 'knex';
import { router } from './association.js';

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

describe('queryAssociation', () => {
  test('should return an error if the required parameters study, cancer, and strategy are not provided', async () => {
    const response = await request(app).get('/signature_association').query({});

    expect(response.statusCode).toBe(400);
    expect(response.body).toEqual({
      error:
        'Missing one or more of the following parameters: study, cancer, strategy',
    });
  });

  test('should query the Signature table with the provided parameters and return the results', async () => {
    const query = {
      study: 'PCAWG',
      strategy: 'WGS',
      cancer: 'Biliary-AdenoCA',
      limit: 10,
    };

    const response = await request(app)
      .get('/signature_association')
      .query(query);

    expect(response.statusCode).toBe(200);
    expect(response.body.length).toBeGreaterThan(0);
  });
});
