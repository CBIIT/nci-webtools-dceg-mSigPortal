export default function TMB(data, study) {
  const genome = { PCAWG: 'GRCh37', TCGA: 'GRCh37' };
  const genomeSize = { GRCh37: 3101976562 / Math.pow(10, 6) };

  const groupByCancer = data.reduce((groups, e) => {
    const cancerName = e.Cancer_Type;
    groups[cancerName] = groups[cancerName] ? [...groups[cancerName], e] : [e];
    return groups;
  }, {});

  const groupBySamples = Object.entries(groupByCancer).map(
    ([cancerName, cancerArray]) => {
      const sampleGroups = cancerArray.reduce((sampleGroup, e) => {
        const sampleName = e.Sample;
        sampleGroup[sampleName] = sampleGroup[sampleName]
          ? [...sampleGroup[sampleName], e]
          : [e];
        return sampleGroup;
      }, {});

      return { cancer: cancerName, samples: sampleGroups };
    },
    {}
  );

  const burdenQuotient = genomeSize[genome[study]];
  const cancerBurden = groupBySamples.map(({ cancer, samples }) => {
    const sampleBurden = Object.entries(samples).map(
      ([sampleName, sampleArray]) => ({
        sample: sampleName,
        burden: Math.log10(
          sampleArray.reduce((sum, e) => e.Exposure + sum, 0) / burdenQuotient
        ),
      })
    );

    return { cancer, sampleBurden };
  });

  console.log(cancerBurden);

  return { traces: [], layout: {} };
}
