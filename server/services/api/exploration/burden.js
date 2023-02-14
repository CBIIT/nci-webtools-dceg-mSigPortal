import { groupBy } from 'lodash-es';

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
// Calculate the number of mutations (burden) per megabase for each study
function addBurden(data, study) {
  // calculate burden per cancer/sample and append to data
  const groupByCancer = groupBy(data, 'cancer');
  const cancerBurden = Object.entries(groupByCancer)
    .map(([cancer, values]) => {
      const groupBySample = groupBy(values, 'sample');
      const samples = Object.entries(groupBySample).reduce(
        (acc, [sample, e]) => {
          acc[sample] = calculateBurden(
            e.reduce((sum, e) => e.exposure + sum, 0),
            study
          );
          return acc;
        },
        {}
      );
      return samples;
    })
    .reduce((obj, e) => ({ ...obj, ...e }), {});
  return data.map((e) => ({
    ...e,
    cancerBurden: cancerBurden[e.sample],
    burden: calculateBurden(e.exposure, study),
  }));
}

export { addBurden };
