import { Worker } from 'worker_threads';
import { fileURLToPath } from 'url';

const workerPath = fileURLToPath(
  new URL('./treeLeafLayout.worker.js', import.meta.url)
);

/**
 * Computes tree/leaf node and link positions in a worker thread so the layout's
 * O(n^2) simulation doesn't block the server's event loop.
 * @param {object} params
 * @param {object} params.data hierarchy ({ name, children })
 * @param {Map<string, number>} params.mutations sample name -> mutation count
 * @param {number} [params.radius]
 * @returns {Promise<{ nodes: object[], links: object[] }>}
 */
export function computeTreeLeafLayout(params) {
  return new Promise((resolve, reject) => {
    const worker = new Worker(workerPath, { workerData: params });
    worker.once('message', ({ layout, error }) => {
      if (error) reject(new Error(error));
      else resolve(layout);
    });
    worker.once('error', reject);
    worker.once('exit', (code) => {
      if (code !== 0) {
        reject(new Error(`Tree/leaf layout worker exited with code ${code}`));
      }
    });
  });
}
