// worker_threads entry point: computes the tree/leaf force layout off the main event loop.
import { parentPort, workerData } from 'worker_threads';
import { createForceDirectedTree } from './treeLeafLayout.js';

try {
  parentPort.postMessage({ layout: createForceDirectedTree(workerData) });
} catch (error) {
  parentPort.postMessage({ error: error.message });
}
