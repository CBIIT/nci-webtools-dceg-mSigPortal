import path from 'path';
import express from 'express';
import Router from 'express-promise-router';
import { check } from 'express-validator';
import multer from 'multer';
import { submit, query, run } from './analysis.js';
import {
  handleValidationErrors,
  logRequests,
  logErrors,
  logFiles,
} from './middleware.js';
import DiskStorage from './storage.js';

export function createApi(env) {
  // define middleware
  const storage = new DiskStorage({
    filename: (req, file) => file.originalname,
    destination: (req) => path.resolve(env.INPUT_FOLDER, req.params.id),
  });
  const upload = multer({ storage });
  const validate = check('id').isUUID();

  // register middleware
  const router = Router();
  router.use(express.json());
  router.use(logRequests());

  // serve static files under /data
  router.use('/data', express.static(env.DATA_FOLDER));

  router.get('/ping', async (req, res) => res.json(true));

  router.post(
    '/upload/:id',
    validate,
    handleValidationErrors,
    upload.any(),
    logFiles(),
    (req, res) => res.json({ id: req.params.id })
  );

  router.post(
    '/submit/:id',
    validate,
    handleValidationErrors,
    async (req, res) =>
      res.json(await submit({ ...req.body, id: req.params.id }, req.app))
  );

  router.get('/run/:id', validate, handleValidationErrors, async (req, res) =>
    res.json(await run(req.params.id, req.app, env))
  );

  router.post('/query/:id', async (req, res) =>
    res.json(await query({ ...req.body, id: req.params.id }))
  );
  router.use(logErrors());
  return router;
}
