export default function TMB(data, study) {
  const genome = { PCAWG: "GRCh37", TCGA: "GRCh37" };
  const genomeSize = { GRCh37: 3101976562 / Math.pow(10, 6) };

  function average(arr) {
    const sum = arr.reduce((a, b) => a + b, 0);
    return sum / arr.length || 0;
  }
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
  // const traces = Object.entries(cancerBurden).map(
  //   ([cancer, samples], groupIndex, array) => ({
  //     cancer: samples.cancer,
  //     samples: samples,
  //     groupIndex: groupIndex,
  //     array: array,
  //     type: "scatter",
  //     marker: { symbol: "circle-open" },
  //     mode: "markers",
  //     testArray: samples.sampleBurden.map((e, i) => array.slice(0, groupIndex)),
  //     testArray2: samples.sampleBurden.map((e, i) =>
  //       array
  //         .slice(0, groupIndex)
  //         .reduce((x0, [_, sigs]) => x0 + sigs.length, 0)
  //     ),
  //     x: samples.sampleBurden.map(
  //       (e, i) => groupIndex * samples.sampleBurden.length + i
  //     ),
  //     y: samples.sampleBurden.map((e) => e.burden),
  //   })
  // );
  const traces = cancerBurden.map((element, index, array) => ({
    element: element,
    index: index,
    array: array,
    name: element.cancer,
    type: "scatter",
    marker: { symbol: "circle-open" },
    mode: "markers",
    // x: element.sampleBurden.map(
    //   (e, i) => index * element.sampleBurden.length + i
    // ),
    y: element.sampleBurden.map((e) => e.burden),
    average: average(element.sampleBurden.map((e) => e.burden)),
    test0: element.sampleBurden.map((e, i) => array.slice(0, index)),
    testx0: element.sampleBurden.map((e, i) =>
      array.slice(0, index).reduce((x0, curr) => x0, 0)
    ),
    testcurr: element.sampleBurden.map((e, i) =>
      array.slice(0, index).reduce((x0, curr) => curr, 0)
    ),
    test1: element.sampleBurden.map((e, i) =>
      array.slice(0, index).reduce((x0, curr) => curr.sampleBurden.length, 0)
    ),
    test2: element.sampleBurden.map((e, i) =>
      array.slice(0, index).reduce((x0, curr, currIndex, a) => x0, 0)
    ),
    x: element.sampleBurden.map(
      (e, i) =>
        array
          .slice(0, index)
          .reduce((x0, curr) => x0 + curr.sampleBurden.length, 0) + i
    ),
  }));
  console.log("traces:--");
  console.log(traces);

  // const shapes = Object.entries(cancerBurden).map(
  //   ([mutation, signatures], groupIndex, array) => ({
  //     type: "rect",
  //     xref: "x",
  //     yref: "paper",
  //     x0: array
  //       .slice(0, groupIndex)
  //       .reduce((x0, [_, signatures]) => x0 + signatures.length, -0.4),
  //     // x0: groupIndex * 16 - 0.4,
  //     x1: array
  //       .slice(0, groupIndex + 1)
  //       .reduce((x0, [_, signatures]) => x0 + signatures.length, -0.6),
  //     // x1: groupIndex * 16 + signatures.length - 0.6,
  //     y0: 1,
  //     y1: 0,

  //     line: {
  //       width: 0,
  //     },
  //     opacity: 0.2,
  //   })
  // );
  const shapes = cancerBurden.map((element, index, array) => ({
    type: "rect",
    xref: "x",
    yref: "paper",
    // x0: array
    //   .slice(0, index)
    //   .reduce((x0, [_, signatures]) => x0 + signatures.length, -0.4),
    x0: index * 38 - 0.4,
    // x1: array
    //   .slice(0, index + 1)
    //   .reduce((x0, [_, signatures]) => x0 + signatures.length, -0.6),
    x1: index * 38 + element.sampleBurden.length - 0.6,
    y0: 0,
    y1: 1,
    fillcolor: "#d3d3d3",
    line: {
      width: 0,
    },
    opacity: 0.5,
  }));

  console.log("shapes:--");
  console.log(shapes);

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

    shapes: shapes,
    //annotations: annotations,
  };
  console.log("layout:");
  console.log(layout);

  return { traces: [...traces], layout: { layout } };
}
