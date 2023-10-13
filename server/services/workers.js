import path from 'path';
import ECS, { ECSClient, RunTaskCommand } from '@aws-sdk/client-ecs';
import Batch, {
  BatchClient,
  CancelJobCommand,
  SubmitJobCommand,
} from '@aws-sdk/client-batch';
import { readJson } from './utils.js';
import { createLogger } from '../services/logger.js';
import { profilerExtraction } from './api/visualization/profilerExtraction.js';

export function getWorker(workerType = 'local') {
  switch (workerType) {
    case 'local':
      return runLocalWorker;
    case 'fargate':
      return runFargateWorker;
    case 'batch':
      return runBatchWorker;
    default:
      throw new Error(`Unknown worker type: ${workerType}`);
  }
}

/**
 * Executes a worker process locally.
 * @param {string} id
 * @param {string} app
 * @param {string} taskName
 * @param {object} env
 * @returns
 */
export async function runLocalWorker(id, app, taskName, env = process.env) {
  const paramsFilePath = path.resolve(env.INPUT_FOLDER, id, 'params.json');
  const params = await readJson(paramsFilePath);
  const logger = app.locals.logger;
  const dbConnection = app.locals.sqlite(id, 'local');

  if (taskName === 'visualization') {
    return await profilerExtraction(params, logger, dbConnection, env);
  } else if (taskName === 'extraction') {
    fetch(`${env.API_BASE_URL}/extraction/run/${id}`);
  }
}

/**
 * Executes a worker process in an AWS Fargate task.
 * @param {string} id
 * @param {string} app
 * @param {string} taskName
 * @param {object} env
 * @returns {Promise<ECS.RunTaskCommandOutput>} task output
 */
export async function runFargateWorker(id, app, taskName, env = process.env) {
  const { ECS_CLUSTER, SUBNET_IDS, SECURITY_GROUP_IDS } = env;
  const taskTypes = {
    visualization: { name: 'worker', taskDefinition: env.WORKER_TASK_NAME },
    extraction: {
      name: 'extraction-worker',
      taskDefinition: env.EXTRACTION_WORKER_TASK_NAME,
    },
  };
  const client = new ECSClient();
  const workerCommand = [
    'node',
    '--require',
    'dotenv/config',
    '--max-old-space-size=16384',
    'worker.js',
    id,
  ];
  const logger = createLogger(env.APP_NAME, env.LOG_LEVEL);
  const taskCommand = new RunTaskCommand({
    cluster: ECS_CLUSTER,
    count: 1,
    launchType: 'FARGATE',
    networkConfiguration: {
      awsvpcConfiguration: {
        securityGroups: SECURITY_GROUP_IDS.split(','),
        subnets: SUBNET_IDS.split(','),
      },
    },
    taskDefinition: taskTypes[taskName].taskDefinition,
    overrides: {
      containerOverrides: [
        {
          name: taskTypes[taskName].name,
          command: workerCommand,
        },
      ],
    },
  });
  const response = await client.send(taskCommand);
  logger.info('Submitted Fargate RunTask command');
  logger.info(workerCommand);
  logger.info(response);
  return response;
}

/**
 * Executes a worker process in an AWS Batch task.
 * @param {string} id
 * @param {string} app
 * @param {string} taskName
 * @param {object} env
 * @returns {Promise<Batch.>} task output
 */
export async function runBatchWorker(id, app, taskName, env = process.env) {
  const { BATCH_JOB_QUEUE, BATCH_JOB_DEFINITION } = env;
  const client = new BatchClient();
  const workerCommand = ['node', '--require', 'dotenv/config', 'worker.js', id];
  const logger = createLogger(env.APP_NAME, env.LOG_LEVEL);
  const jobCommand = new SubmitJobCommand({
    // SubmitJobRequest
    jobName: `${taskName}-${id}`,
    jobQueue: BATCH_JOB_QUEUE,
    jobDefinition: BATCH_JOB_DEFINITION,
    containerOverrides: {
      command: workerCommand,
    },
    propagateTags: true,
  });
  const response = await client.send(jobCommand);
  logger.info('Submitted Batch SubmitJob command');
  logger.info(workerCommand);
  logger.info(response);
  return response;
}
