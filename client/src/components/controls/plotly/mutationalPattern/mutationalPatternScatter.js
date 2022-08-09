export default function mutationalPatternScatter(
  data,
  type,
  pattern1,
  pattern2
) {
  const maxTotal = Math.max(...data.map((o) => o.total));
  console.log(maxTotal);

  var trace1 = {
    x: data.map((e) => e.n1),
    y: data.map((e) => e.n2),
    customdata: data.map((e) => ({ total: e.total })),
    mode: 'markers',
    type: 'scatter',
    opacity: 1,
    marker: {
      color: 'white',
      line: {
        color: 'black',
        width: 1,
      },
      size: data.map((e, i, a) =>
        e.total < 1000
          ? 6
          : e.total < 2000
          ? 7
          : e.total < 3000
          ? 8
          : e.total < 4000
          ? 9
          : e.total < 5000
          ? 10
          : e.total < 10000
          ? 13
          : e.total < 100000
          ? 15
          : e.total < 1000000
          ? 17
          : 20
      ),
    },
    hovertemplate:
      '<b>' +
      pattern2 +
      ':</b>' +
      ' %{x} <br>' +
      '<b>' +
      pattern1 +
      ':</b>' +
      ' %{y}<br>' +
      '<b>Total:</b> %{customdata.total}<extra></extra>',
    showlegend: true,
    legendgroup: 'size',
    legendgrouptitle: {
      text: 'Number of mutations',
    },
  };
  var trace2 = {
    name: data[0].study + '@' + data[0].cancer,
    x: data.map((e) => e.n1),
    y: data.map((e) => e.n2),
    customdata: data.map((e) => ({ total: e.total })),
    mode: 'markers',
    type: 'scatter',
    opacity: 1,
    marker: {
      color: 'green',
      line: {
        color: 'black',
        width: 1,
      },
      size: data.map((e, i, a) =>
        e.total < 1000
          ? 6
          : e.total < 2000
          ? 7
          : e.total < 3000
          ? 8
          : e.total < 4000
          ? 9
          : e.total < 5000
          ? 10
          : e.total < 10000
          ? 13
          : e.total < 100000
          ? 15
          : e.total < 1000000
          ? 17
          : 20
      ),
    },
    hovertemplate:
      '<b>' +
      pattern2 +
      ':</b>' +
      ' %{x} <br>' +
      '<b>' +
      pattern1 +
      ':</b>' +
      ' %{y}<br>' +
      '<b>Total:</b> %{customdata.total}<extra></extra>',
    showlegend: true,
    legendgrouptitle: {
      text: 'Study',
    },
  };

  console.log(trace1);

  console.log(trace2);

  const traces = [trace1, trace2];
  var layout = {
    height: 700,
    hoverlabel: { bgcolor: '#FFF' },
    legend: {},
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
