import request from 'supertest';
import express from 'express';
import knex from 'knex';
import { router as signatureRouter } from './signature.js';

const env = process.env;
const app = express();

app.use(express.json());
app.use(signatureRouter);

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

describe('querySignature', () => {
  test('should return an error if the required parameters study, cancer, and strategy are not provided', async () => {
    const response = await request(app).get('/mutational_signature').query({});

    expect(response.statusCode).toBe(400);
    expect(response.body).toEqual({
      error:
        'Missing one or more of the following parameters: study, cancer, strategy',
    });
  });

  test('should query the Signature table with the provided parameters and return the results', async () => {
    const source = 'Reference_signatures';
    const strategy = 'WGS';
    const profile = 'SBS';
    const matrix = '96';
    const limit = 10;

    const response = await request(app)
      .get('/mutational_signature')
      .query({ source, strategy, profile, matrix, limit });

    expect(response.statusCode).toBe(200);
    expect(response.body.length).toBeGreaterThan(0);
  });
});

describe('signatureOptions', () => {
  test('should return an error if the required parameters userId or study are not provided', async () => {
    const response = await request(app)
      .get('/mutational_signature_options')
      .query({});

    expect(response.statusCode).toBe(400);
    expect(response.body).toEqual({
      error: 'Missing one or more of the following parameters: userId or study',
    });
  });

  test('should query the signature_options table with the provided parameters and return the results', async () => {
    const source = 'Reference_signatures';
    const strategy = 'WGS';
    const profile = 'SBS';
    const limit = 10;

    const response = await request(app)
      .get('/mutational_signature_options')
      .query({ source, strategy, profile, limit });

    expect(response.statusCode).toBe(200);
    expect(response.body.length).toBeGreaterThan(0);
  });
});

describe('signatureSummary', () => {
  test('should query the signature_summary table with the provided parameters and return the results', async () => {
    const profile = 'SBS';
    const matrix = '96';
    const limit = 10;

    const response = await request(app)
      .get('/mutational_signature_summary')
      .query({ profile, matrix, limit });

    expect(response.statusCode).toBe(200);
    expect(response.body.length).toBeGreaterThan(0);
  });
});
