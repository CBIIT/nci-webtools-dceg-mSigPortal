import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';
import { after, before, test } from 'node:test';
import { createServer } from 'vite';

let server;
let utils;
let MsAssociation;

before(async () => {
  // Load the actual plot through Vite so JSX and the app's alias work without
  // duplicating its implementation or adding another test dependency.
  server = await createServer({
    configFile: false,
    root: fileURLToPath(new URL('..', import.meta.url)),
    resolve: {
      alias: { '@': fileURLToPath(new URL('../src', import.meta.url)) },
    },
    server: { middlewareMode: true, watch: null, hmr: false, ws: false },
    optimizeDeps: { noDiscovery: true, include: [] },
    appType: 'custom',
  });
  utils = await server.ssrLoadModule('/src/components/controls/utils/utils.js');
  ({ default: MsAssociation } = await server.ssrLoadModule(
    '/src/components/controls/plotly/msAssociation/msAssociation.jsx'
  ));
});

after(async () => server?.close());

function closeTo(actual, expected, tolerance = 1e-10) {
  assert.ok(
    Number.isFinite(actual) && Math.abs(actual - expected) <= tolerance,
    `Expected ${actual} to be within ${tolerance} of ${expected}`
  );
}

function exposureRows(x, y) {
  return x.flatMap((exposure, i) => [
    { sample: `sample-${i}`, signatureName: 'SBS7a', exposure },
    { sample: `sample-${i}`, signatureName: 'SBS2', exposure: y[i] },
  ]);
}

function annotation(plot) {
  return plot.layout.annotations[0].text;
}

test('round supports scientific notation and preserves decimal rounding', () => {
  assert.equal(utils.round(1.6525311544709337e-19, 3), 0);
  assert.equal(utils.round(1.65e-19, 21), 1.65e-19);
  assert.equal(utils.round(1.005, 2), 1.01);
  assert.equal(utils.round(-1.1480147617436247, 2), -1.15);
  assert.equal(utils.round(1e25, 2), 1e25);
  assert.equal(utils.round(1e308, 2), 1e308);
  assert.equal(utils.round(Infinity, 2), Infinity);
});

test('Pearson statistic, two-sided p, and Fisher-z CI match R cor.test', () => {
  const stats = utils.calculatePearson([1, 2, 3, 4, 5], [2, 1, 4, 3, 6]);
  closeTo(stats.pcorr, 0.8219949365267864);
  closeTo(stats.statistic, 2.5);
  closeTo(stats.pValue, 0.08770664700806559);
  closeTo(stats.CILower, -0.2193433574707815);
  closeTo(stats.CIUpper, 0.987852453504985);
});

test('Spearman handles tied ranks and matches R asymptotic p', () => {
  const stats = utils.calculateSpearman([0, 0, 0, 1, 2, 2], [1, 0, 0, 2, 2, 3]);
  closeTo(stats.rho, 0.9057894597833128);
  closeTo(stats.pValue, 0.01289535007704897);
  // R cor.test(rank(x), rank(y)) provides the approximate Fisher-z interval.
  closeTo(stats.CILower, 0.3557171776825835);
  closeTo(stats.CIUpper, 0.9897683843103593);
});

test('Sherlock screenshot regression: correct Spearman rho and p', async () => {
  const text = await readFile(
    new URL(
      '../public/assets/exampleInput/Sherlock_SBS96_exposure.txt',
      import.meta.url
    ),
    'utf8'
  );
  const [header, ...rows] = text
    .trim()
    .split(/\r?\n/)
    .map((row) => row.split('\t'));
  const x = rows.map((row) => Number(row[header.indexOf('SBS7a')]));
  const y = rows.map((row) => Number(row[header.indexOf('SBS2')]));
  const stats = utils.calculateSpearman(
    x.map((v) => Math.log10(v + 1)),
    y.map((v) => Math.log10(v + 1))
  );
  assert.equal(stats.n, 260);
  closeTo(stats.rho, -0.066023545147264);
  closeTo(stats.pValue, 0.2888599995130595);
  const textAnnotation = annotation(
    MsAssociation(exposureRows(x, y), 'SBS7a', 'SBS2')
  );
  assert.match(textAnnotation, /Pearson:.*p = 0\.252/);
  assert.match(textAnnotation, /Spearman:.*p \(approx\.\) = 0\.289/);
  assert.doesNotMatch(textAnnotation, /NaN|Infinity/);
});

test('constant inputs return unavailable statistics for both methods', () => {
  for (const constant of [0, 0.1]) {
    for (const [x, y] of [
      [Array(38).fill(constant), Array.from({ length: 38 }, (_, i) => i)],
      [Array.from({ length: 38 }, (_, i) => i), Array(38).fill(constant)],
      [Array(38).fill(constant), Array(38).fill(constant)],
    ]) {
      const pearson = utils.calculatePearson(x, y);
      const spearman = utils.calculateSpearman(x, y);
      for (const stats of [pearson, spearman]) {
        assert.match(stats.reason, /no variation/);
        assert.equal(stats.pValue, null);
        assert.equal(stats.CILower, null);
        assert.equal(stats.CIUpper, null);
      }
      assert.equal(pearson.pcorr, null);
      assert.equal(spearman.rho, null);
    }
  }
});

test('PCAWG screenshot pattern: 38 zero SBS7a exposures show an explanation', () => {
  const data = exposureRows(
    Array(38).fill(0),
    Array.from({ length: 38 }, (_, i) => i)
  );
  const plot = MsAssociation(data, 'SBS7a', 'SBS2');
  const text = annotation(plot);
  assert.match(
    text,
    /Pearson: N\/A \(one or both signatures have no variation\)/
  );
  assert.match(
    text,
    /Spearman: N\/A \(one or both signatures have no variation\)/
  );
  assert.match(text, /n<sub>pairs<\/sub> = 38/);
  assert.doesNotMatch(text, /NaN|Infinity/);
  assert.ok(plot.traces.every((trace) => trace.mode !== 'lines'));
});

test('very small nonzero p values stay in scientific notation', () => {
  const x = Array.from({ length: 30 }, (_, i) => i + 1);
  const y = x.map((v, i) => v + (i % 2 ? 0.75 : -0.75));
  const text = annotation(MsAssociation(exposureRows(x, y), 'SBS7a', 'SBS2'));
  assert.match(text, /p = \d\.\d{2}e-\d+/);
  assert.match(text, /p \(approx\.\) = \d\.\d{2}e-\d+/);
  assert.doesNotMatch(text, /NaN|Infinity|p = 0[, ]/);
});

test('perfect positive and negative correlation do not produce NaN', () => {
  for (const y of [
    [0, 1, 2, 3],
    [3, 2, 1, 0],
  ]) {
    for (const calculate of [utils.calculatePearson, utils.calculateSpearman]) {
      const stats = calculate([0, 1, 2, 3], y);
      const r = stats.pcorr ?? stats.rho;
      assert.equal(Math.abs(r), 1);
      assert.equal(stats.pValue, 0);
      assert.equal(stats.CILower, r);
      assert.equal(stats.CIUpper, r);
    }
  }
  const plot = MsAssociation(
    exposureRows([0, 1, 3, 7], [0, 1, 3, 7]),
    'SBS7a',
    'SBS2'
  );
  assert.match(annotation(plot), /t<sub>Student<\/sub> = ∞/);
  assert.match(annotation(plot), /p = &lt; 0\.001/);
  assert.doesNotMatch(annotation(plot), /NaN|Infinity/);
});

test('empty, single-sample, and both-signatures-filtered plots explain missing statistics', () => {
  const plots = [
    MsAssociation([], 'SBS7a', 'SBS2'),
    MsAssociation(exposureRows([1], [2]), 'SBS7a', 'SBS2'),
    MsAssociation(exposureRows([0, 0], [1, 2]), 'SBS7a', 'SBS2', true),
  ];
  for (const plot of plots) {
    assert.match(annotation(plot), /at least two paired samples are required/);
    assert.doesNotMatch(annotation(plot), /NaN|Infinity/);
  }
});

test('two pairs have no test or CI; three pairs have a test but no CI', () => {
  for (const calculate of [utils.calculatePearson, utils.calculateSpearman]) {
    const two = calculate([0, 1], [1, 2]);
    assert.equal(two.pValue, null);
    assert.equal(two.CILower, null);
    assert.match(two.testReason, /three/);
    assert.match(two.ciReason, /four/);
    const three = calculate([0, 1, 2], [1, 0, 2]);
    closeTo(three.pValue, 2 / 3, 1e-8);
    assert.equal(three.CILower, null);
    assert.match(three.ciReason, /four/);
  }
  for (const data of [
    exposureRows([0, 1], [1, 2]),
    exposureRows([0, 1, 2], [1, 0, 2]),
  ]) {
    assert.doesNotMatch(
      annotation(MsAssociation(data, 'SBS7a', 'SBS2')),
      /NaN|Infinity/
    );
  }
});

test('invalid and mismatched correlation inputs fail with clear errors', () => {
  for (const calculate of [utils.calculatePearson, utils.calculateSpearman]) {
    assert.throws(() => calculate([1, 2], [1]), /same number/);
    assert.throws(() => calculate([1, NaN], [1, 2]), /finite numeric/);
    assert.throws(() => calculate([1, 2], [1, Infinity]), /finite numeric/);
  }
});
