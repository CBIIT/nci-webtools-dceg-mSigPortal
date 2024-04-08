// import request from 'supertest';
// import express from 'express';
// import knex from 'knex';
// import { router } from './visualization/visualization.js';

// const env = process.env;
// const app = express();
// app.use(express.json());
// app.use('/', router);

// app.locals.connection = knex({
//   client: 'postgres',
//   connection: {
//     host: env.POSTGRES_HOST,
//     port: env.POSTGRES_PORT,
//     user: env.POSTGRES_USER,
//     password: env.POSTGRES_PASS,
//     database: env.POSTGRES_DB,
//   },
// });

// describe('Visualization API', () => {
//   test('GET /mutational_spectrum', async () => {
//     const response = await request(app).get('/mutational_spectrum');
//     expect(response.statusCode).toBe(200);
//     // Add more assertions based on your expected response
//   });

//   test('GET /mutational_spectrum_options', async () => {
//     const response = await request(app).get('/mutational_spectrum_options');
//     expect(response.statusCode).toBe(200);
//     // Add more assertions based on your expected response
//   });

//   test('GET /mutational_spectrum_summary', async () => {
//     await request(app)
//       .get('/mutational_spectrum_summary')
//       .query({ study: 'PCAWG', cancer: 'Lung-AdenoCA', strategy: 'WGS' })
//       .expect(200);
//     // Add more assertions based on your expected response
//   });

//   // user data only
//   // test('GET /cluster', async () => {
//   //   const response = await request(app).get('/cluster');
//   //   expect(response.statusCode).toBe(200);
//   //   // Add more assertions based on your expected response
//   // });

//   test('GET /refgenome', async () => {
//     const response = await request(app).get('/refgenome');
//     expect(response.statusCode).toBe(200);
//     // Add more assertions based on your expected response
//   });
// });
