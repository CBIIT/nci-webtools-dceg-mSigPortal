export default function TMB(data) {
  console.log(data);
  function average(arr) {
    const sum = arr.reduce((a, b) => a + b, 0);
    return sum / arr.length || 0;
  }

  const totalCancer = data.length;
  console.log(totalCancer);

  const absYValue = data
    .map((o) => o.samples.map((e) => Math.abs(e.burden)))
    .flat();
  const yMax = Math.max(...absYValue);

  const traces = data.map((element, index, array) => ({
    // element: element,
    // index: index,
    // array: array,
    // name: `${element.cancer}`,
    type: 'scatter',
    marker: { symbol: 'circle-open', size: 3, color: 'black' },
    mode: 'markers',
    y: element.samples.map((e) => e.burden),
    // average: average(element.samples.map((e) => e.tmb)),
    hovertemplate: 'Number of mutations: %{y}<br>',
    // x: element.samples.map(
    //   (e, i) =>
    //     array
    //       .slice(0, index)
    //       .reduce((x0, curr) => x0 + curr.samples.length, 0) +
    //     i +
    //     0.5
    // ),

    x: element.samples.map(
      (e, i) => index + 0.1 + (0.8 / element.samples.length) * i
    ),
  }));

  console.log('traces:--');
  console.log(traces);

  const topLabel = data.map((element, index, array) => ({
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    // x:
    //   array
    //     .slice(0, index)
    //     .reduce((x0, curr) => x0 + curr.samples.length, 0) +
    //   (element.samples.length - 1) * 0.5,
    //x: array.length * 0.5,
    x: array.length > 1 ? index : (index + index + 1) * 0.5,
    y: 1.01,
    text: `${element.cancer}`,
    showarrow: false,
    font: {
      size: 10,
    },
    align: 'center',
    textangle: 55,
  }));
  console.log('top label:--');
  console.log(topLabel);

  const bottoLabel1 = data.map((element, index, array) => ({
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: (index + index + 1) * 0.5,
    y: -0.1,
    text: `${element.samples.length}`,
    showarrow: false,
    font: {
      size: 12,
      color: 'blue',
    },
    align: 'center',
  }));

  const bottoLabelline = data.map((element, index, array) => ({
    type: 'line',
    xref: 'x',
    yref: 'paper',
    x0: index + 0.4,
    x1: index + 0.6,
    y0: -0.11,
    y1: -0.11,
    line: {
      width: 1,
      color: 'black',
    },
  }));

  const bottoLabel2 = data.map((element, index, array) => ({
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: (index + index + 1) * 0.5,
    y: -0.18,
    text: `${element.samples.length}`,
    showarrow: false,
    font: {
      size: 12,
      color: 'green',
    },
    align: 'center',
  }));

  const shapes = data.map((element, index, array) => ({
    type: 'rect',
    xref: 'x',
    yref: 'paper',

    // x0: array
    //   .slice(0, index)
    //   .reduce((x0, curr) => x0 + curr.samples.length, 0),

    // x1: array
    //   .slice(0, index + 1)
    //   .reduce((x0, curr) => x0 + curr.samples.length, 0),
    x0: index,
    x1: index + 1,
    y0: 0,
    y1: 1,
    fillcolor: index % 2 === 0 ? 'gray' : '#F8F8F8',
    line: {
      width: 0,
    },
    opacity: 0.2,
  }));
  console.log('shapes:--');
  console.log(shapes);

  const lines = data.map((element, index, array) => ({
    type: 'line',
    xref: 'x',
    yref: 'y',

    // x0: array
    //   .slice(0, index)
    //   .reduce((x0, curr) => x0 + curr.samples.length, 0),

    // x1: array
    //   .slice(0, index + 1)
    //   .reduce((x0, curr) => x0 + curr.samples.length, 0),
    x0: index + 0.1,
    x1: index + 0.9,

    y0: element.medianBurden,
    y1: element.medianBurden,
    line: {
      width: 1,
      color: 'red',
    },
  }));

  const layout = {
    // title: {
    //   text: "Tumor Mutational Burden Separated by Signatures",
    //   yanchor: "top",
    // },
    width: totalCancer > 1 ? null : 350,
    autosize: true,
    height: 500,
    showlegend: false,
    xaxis: {
      showticklabels: false,
      tickfont: {
        size: 10,
      },
      autorange: false,
      range: [0, totalCancer],
      linecolor: 'black',
      linewidth: 2,
      mirror: true,
      tickmode: 'array',
      showgrid: false,
      //tickvals: flatSorted.map((_, i) => i),
      //ticktext: flatSorted.map((_, i) => i),
    },
    yaxis: {
      title: 'Number of Mutations per Megabase<br>(log10)',
      zeroline: false,
      //showline: true,
      linecolor: 'black',
      linewidth: 2,
      mirror: true,
      automargin: true,
      autorange: true,
      range: [-Math.floor(yMax), Math.floor(yMax)],
    },

    shapes: [...shapes, ...lines, ...bottoLabelline],
    annotations: [...topLabel, ...bottoLabel1, ...bottoLabel2],
  };

  console.log('layout:');
  console.log(layout);

  var config = {
    //responsive: true,
  };

  return { traces: [...traces], layout: layout, config };
  //return { traces, layout };
}
