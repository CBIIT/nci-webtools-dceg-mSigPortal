import request from 'supertest';
import express from 'express';
import knex from 'knex';
import { router as visualizationRouter } from './visualization.js';

const env = process.env;
const app = express();

app.use(express.json());
app.use(visualizationRouter);

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

describe('Seqmatrix API', () => {
  // beforeAll(async () => {});

  // afterAll(async () => {
  //   await connection.destroy();
  // });

  test('Should return an error if the required parameters study, cancer, and strategy are not provided', async () => {
    const response = await request(app).get('/mutational_spectrum').query({});

    expect(response.statusCode).toBe(400);
    expect(response.body).toEqual({
      error:
        'Missing one or more of the following parameters: study, cancer, strategy',
    });
  });

  test('Should query the Seqmatrix table with the provided parameters and return the results', async () => {
    const study = 'Sherlock-Lung-232';
    const cancer = 'LCINS';
    const strategy = 'WGS';
    const limit = 10;

    const response = await request(app)
      .get('/mutational_spectrum')
      .query({ study, cancer, strategy, limit });

    expect(response.statusCode).toBe(200);
    expect(response.body.length).toBeGreaterThan(0);
  });
});

describe('seqmatrixOptions', () => {
  test('should return an error if the required parameters study, cancer, and strategy are not provided', async () => {
    const response = await request(app)
      .get('/mutational_spectrum_options')
      .query({});

    expect(response.statusCode).toBe(400);
    expect(response.body).toEqual({
      error:
        'Missing one or more of the following parameters: study, cancer, strategy',
    });
  });

  test('should query the seqmatrix_options table with the provided parameters and return the results', async () => {
    const study = 'Sherlock-Lung-232';
    const cancer = 'LCINS';
    const strategy = 'WGS';
    const limit = 10;

    const response = await request(app)
      .get('/mutational_spectrum_options')
      .query({ study, cancer, strategy, limit });

    expect(response.statusCode).toBe(200);
    expect(response.body.length).toBeGreaterThan(0);
  });
});

describe('seqmatrixSummary', () => {
  test('should return an error if the required parameters study, cancer, and strategy are not provided', async () => {
    const response = await request(app)
      .get('/mutational_spectrum_summary')
      .query({});

    expect(response.statusCode).toBe(400);
    expect(response.body).toEqual({
      error:
        'Missing one or more of the following parameters: study, cancer, strategy',
    });
  });

  test('should query the seqmatrix_summary table with the provided parameters and return the results', async () => {
    const study = 'Sherlock-Lung-232';
    const cancer = 'LCINS';
    const strategy = 'WGS';
    const limit = 10;

    const response = await request(app)
      .get('/mutational_spectrum_summary')
      .query({ study, cancer, strategy, limit });

    expect(response.statusCode).toBe(200);
    expect(response.body.length).toBeGreaterThan(0);
  });
});

describe('RefGenome API', () => {
  test('Should query the refgenome table with the provided parameters and return the results', async () => {
    const response = await request(app).get('/refgenome').query({});

    expect(response.statusCode).toBe(200);
    expect(response.body.length).toBeGreaterThan(0);
  });
});
