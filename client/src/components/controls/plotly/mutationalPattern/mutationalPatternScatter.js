import { text } from 'd3';

export default function mutationalPatternScatter(
  data,
  type,
  pattern1,
  pattern2
) {
  const maxTotal = Math.max(...data.map((o) => o.total));
  console.log(maxTotal);
  console.log(maxTotal.toExponential());
  function format_output(output) {
    var n = Math.log(output) / Math.LN10;
    var x = Math.round(n);
    if (x < 0) x = 0;
    let out = '1';
    for (var i = 0; i <= x; i++) out += '0';
    return out / 10;
  }
  const maxMutationFilter = format_output(maxTotal);
  console.log(maxMutationFilter);
  const data1 = data.filter((o) => o.total < maxMutationFilter / 100);
  const data2 = data.filter(
    (o) => o.total < maxMutationFilter / 10 && o.total > maxMutationFilter / 100
  );
  const data3 = data.filter((o) => o.total > maxMutationFilter / 10);

  var trace1 = {
    name: (maxMutationFilter / 100).toLocaleString(undefined),
    x: data1.map((e) => e.n1),
    y: data1.map((e) => e.n2),
    mode: 'markers',
    type: 'scatter',
    opacity: 1,
    marker: {
      color: 'white',
      line: {
        color: 'black',
        width: 1,
      },
      size: 5,
    },
    showlegend: true,
    legendgroup: 'size',
    legendgrouptitle: {
      text: 'Number of mutations',
    },
  };
  var trace2 = {
    name: (maxMutationFilter / 10).toLocaleString(undefined),
    x: data2.map((e) => e.n1),
    y: data2.map((e) => e.n2),
    mode: 'markers',
    type: 'scatter',
    opacity: 1,
    marker: {
      color: 'white',
      line: {
        color: 'black',
        width: 1,
      },
      size: 10,
    },
    showlegend: true,
    legendgroup: 'size',
    legendgrouptitle: {
      text: 'Number of mutations',
    },
  };
  var trace3 = {
    name: maxMutationFilter.toLocaleString(undefined),
    x: data3.map((e) => e.n1),
    y: data3.map((e) => e.n2),
    mode: 'markers',
    type: 'scatter',
    opacity: 1,
    marker: {
      color: 'white',
      line: {
        color: 'black',
        width: 1,
      },
      size: 15,
    },
    showlegend: true,
    legendgroup: 'size',
    legendgrouptitle: {
      text: 'Number of mutations',
    },
  };
  var trace4 = {
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
        e.total < maxMutationFilter / 100
          ? 5
          : e.total < maxMutationFilter / 10
          ? 10
          : 15
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

  console.log(trace4);

  const traces = [trace1, trace2, trace3, trace4];
  var layout = {
    height: 700,
    hoverlabel: { bgcolor: '#FFF' },
    legend: {
      x: 1,
      y: 0.5,
    },
    xaxis: {
      title: pattern2,
      range: [-0.1, 1.1],
      zeroline: false,
    },
    yaxis: {
      title: pattern1,
      range: [-0.1, 1.1],
      zeroline: false,
    },
    title: {
      text: '<b>Proportion of Mutational Pattern Context Compared to Other Contexts with the same SBS Mutation</b>',
      font: {
        family: 'Times New Roman',
        size: 18,
      },
    },
  };

  var config = {
    responsive: true,
  };

  return { traces: traces, layout: layout, config };
}
