export default function TMB(data, study) {
  const genome = { PCAWG: "GRCh37", TCGA: "GRCh37" };
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

  console.log(groupByCancer);

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

  const flatSorted = Object.values(cancerBurden).flat();

  console.log(cancerBurden);
  console.log(flatSorted);
  const traces = Object.entries(cancerBurden).map(
    ([cancer, samples], groupIndex, array) => ({
      cancer: samples.cancer,
      samples: samples,
      groupIndex: groupIndex,
      array: array,
      type: "scatter",
      marker: { symbol: "circle-open" },
      mode: "markers",
      testArray: samples.sampleBurden.map((e, i) => array.slice(0, groupIndex)),
      testArray2: samples.sampleBurden
        .map((e, i) => array[i][1].slice(0, groupIndex))
        .reduce((x0, [_, sigs]) => x0 + sigs.length, 0),
      x: samples.sampleBurden.map(
        (e, i) => groupIndex * samples.sampleBurden.length + i
      ),
      y: samples.sampleBurden.map((e) => e.burden),
    })
  );
  console.log("traces:--");
  console.log(traces);
  const annotations = Object.entries(cancerBurden).map(
    ([cancer, samples], groupIndex, array) => ({
      xref: "x",
      yref: "paper",
      x:
        array
          .slice(0, groupIndex)
          .reduce((x0, [_, sigs]) => x0 + sigs.length, 0) +
        (samples.length - 1) * 0.5,
      xanchor: "bottom",
      y: 1.05,
      yanchor: "bottom",
      text: "<b>" + samples.cancer + "</b>",
      showarrow: false,
      font: {
        size: 18,
      },
      align: "center",
    })
  );
  console.log("annotations: --");
  console.log(annotations);

  const layout = {
    xaxis: {
      showticklabels: false,
      showline: true,
      tickangle: -90,
      tickfont: {
        size: 10,
      },
      tickmode: "array",
      //tickvals: flatSorted.map((_, i) => i),
      //ticktext: flatSorted.map((e) => e.mutationType),
    },
    yaxis: {
      title: "Number of Mutations per Megabase",
      autorange: true,
    },

    //shapes: shapes,
    annotations: annotations,
  };
  console.log("layout:");
  console.log(layout);

  return { traces: [...traces], layout: { layout } };
}
