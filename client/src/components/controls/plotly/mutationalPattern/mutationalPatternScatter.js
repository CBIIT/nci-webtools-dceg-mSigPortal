export default function mutationalPatternScatter(
  data,
  type,
  subtype1,
  subtype2,
  pattern1,
  pattern2
) {
  console.log(data);
  console.log(type);

  var scattertrace = {
    x: data.map((e, i, a) => e.n1),
    y: data.map((e, i, a) => e.n2),
    mode: 'markers',
    type: 'scatter',

    marker: { size: 12 },
  };

  console.log(scattertrace);

  const traces = [scattertrace];
  var layout = {
    height: 600,
    xaxis: {
      title: pattern2,
      range: [0, 1],
    },
    yaxis: {
      title: pattern1,
      range: [0, 1],
    },
    title:
      'Proportion of Mutational Pattern Context Compared to Other Contexts with the same SBS Mutation',
  };

  var config = {
    responsive: true,
  };

  return { traces: traces, layout: layout, config };
}
