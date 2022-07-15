const express = require('express');
const logger = require('../services/logger');
const config = require('../config.json');

const {
  visualizationProfilerExtraction,
  getResults,
  visualizationWrapper,
  getSignaturesUser,
  upload,
  visualizationDownload,
  visualizationDownloadPublic,
  explorationWrapper,
  submitQueue,
  getQueueResults,
  getVisExample,
  getExposureExample,
  getPublications,
  getImageS3Batch,
  getImageS3,
  getFileS3,
  downloadWorkspace,
  associationWrapper,
} = require('../services/apiAnalysis');

const visualizationRoutes = require('../services/analysis/visualization/visualization');

const explorationRoutes = require('../services/analysis/exploration/exploration');

const router = express.Router();

router.use(
  '/results',
  express.static(config.results.folder, {
    setHeaders: (res, path, stat) => {
      res.set('Cache-Control', 'max-age=0, must-revalidate');
    },
  })
);

router.post('/profilerExtraction', visualizationProfilerExtraction);

router.post('/getResults', getResults);

router.post('/visualizationWrapper', visualizationWrapper);

router.post('/getSignaturesUser', getSignaturesUser);

router.post('/upload', upload);

router.get('/visualization/download', visualizationDownload);

router.post(
  '/visualization/downloadPublic',

  visualizationDownloadPublic
);

router.post('/explorationWrapper', explorationWrapper);

router.post('/queue', submitQueue);

router.get('/getQueueResults/:id', getQueueResults);

router.get('/getVisExample/:example', getVisExample);

router.get('/getExposureExample/:example', getExposureExample);

router.get('/getPublications', getPublications);

router.post('/getImageS3Batch', getImageS3Batch);

router.post('/getImageS3', getImageS3);

router.post('/getFileS3', getFileS3);

router.post('/downloadWorkspace', downloadWorkspace);

router.post('/associationWrapper', associationWrapper);

// visualization
router.use(visualizationRoutes);

// exploration
router.use(explorationRoutes);

router.use((error, req, res, next) => {
  logger.error(error);
  res.status(500).json(error);
});

module.exports = router;
