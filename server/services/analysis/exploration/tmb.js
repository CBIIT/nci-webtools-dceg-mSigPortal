const { groupBy } = require('lodash');

function calculateBurden(exposure, study) {
  const genome = { PCAWG: 'GRCh37', TCGA: 'GRCh37' };
  const genomeSizeMap = {
    GRCh37: 3101976562 / Math.pow(10, 6),
    GRCh38: 3217346917 / Math.pow(10, 6),
    hg18: 3080436051 / Math.pow(10, 6),
    hg19: 3101976562 / Math.pow(10, 6),
    mmc9: 2654911517 / Math.pow(10, 6),
    mmc10: 2725537669 / Math.pow(10, 6),
  };
  const genomeSize = genomeSizeMap.GRCh37;
  // const genomeSize = genomeSize[genome[study]];
  return Math.log10(exposure / genomeSize);
}

// Calculate the number of mutations per megabase for each study
function calculateTmb(data, study) {
  const groupByCancer = groupBy(data, 'cancer');

  // calculate burden per cancer/sample
  const cancerBurden = Object.entries(groupByCancer)
    .map(([cancer, values]) => {
      const groupBySample = groupBy(values, 'sample');
      const samples = Object.entries(groupBySample).map(([sample, e]) => ({
        sample: sample,
        tmb: calculateBurden(e.reduce((sum, e) => e.exposure + sum, 0)),
      }));
      samples.sort((a, b) => a.tmb - b.tmb);

      // find median tmb
      const tmbs = samples.map((e) => e.tmb);
      const medianTmb =
        tmbs.length % 2 == 0
          ? (tmbs[tmbs.length / 2] + tmbs[tmbs.length / 2 - 1]) / 2
          : tmbs[Math.floor(tmbs.length / 2)];

      return { cancer, samples, medianTmb };
    })
    .sort((a, b) => a.medianTmb - b.medianTmb);

  return cancerBurden;
}

// Calculate the number of mutations per megabase for each signature
function calculateTmbSignature(data, study) {
  const groupBySignature = groupBy(data, 'signatureName');

  // calculate burden per signature/sample
  const signatureBurden = Object.entries(groupBySignature)
    .map(([signatureName, values]) => {
      const groupBySample = groupBy(values, 'sample');
      const samples = Object.entries(groupBySample).map(([sample, e]) => ({
        sample: sample,
        tmb: calculateBurden(e.reduce((sum, e) => e.exposure + sum, 0)),
      }));
      samples.sort((a, b) => a.tmb - b.tmb);

      // find median tmb
      const tmbs = samples.map((e) => e.tmb);
      const medianTmb =
        tmbs.length % 2 == 0
          ? (tmbs[tmbs.length / 2] + tmbs[tmbs.length / 2 - 1]) / 2
          : tmbs[Math.floor(tmbs.length / 2)];

      return { signatureName, samples, medianTmb };
    })
    .sort((a, b) => a.medianTmb - b.medianTmb);

  return signatureBurden;
}

module.exports = { calculateTmb, calculateTmbSignature };
