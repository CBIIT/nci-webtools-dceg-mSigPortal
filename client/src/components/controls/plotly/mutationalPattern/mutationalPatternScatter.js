export default function mutationalPatternScatter(
  data,
  type,
  pattern1,
  pattern2
) {
  const maxSize = 40;
  const maxTotal = Math.max(...data.map((o) => o.total));
  console.log(maxTotal);

  var scattertrace = {
    x: data.map((e, i, a) => e.n1),
    y: data.map((e, i, a) => e.n2),
    mode: 'markers',
    type: 'scatter',

    marker: {
      color: '#499855',
      line: {
        color: 'black',
        width: 0.5,
      },
      size: data.map((e, i, a) =>
        e.total < 1000 ? 4 : e.total < 10000 ? 8 : 12
      ),
    },
  };

  console.log(scattertrace);

  const traces = [scattertrace];
  var layout = {
    height: 700,
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
