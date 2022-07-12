import { groupBy } from 'lodash';

export default function TMB(data, study) {
  // Calculate the number of mutations per megabase for each study
  const genome = { PCAWG: 'GRCh37', TCGA: 'GRCh37' };
  const genomeSize = { GRCh37: 3101976562 / Math.pow(10, 6) };
  const burden = (exposure) => Math.log10(exposure / genomeSize[genome[study]]);

  console.log(data);
  function average(arr) {
    const sum = arr.reduce((a, b) => a + b, 0);
    return sum / arr.length || 0;
  }
  const groupByCancer = groupBy(data, 'cancer');

  // calculate burden per cancer/sample
  const cancerBurden = Object.entries(groupByCancer)
    .map(([cancer, values]) => {
      const groupBySample = groupBy(values, 'sample');
      const tmbs = Object.entries(groupBySample).map(([_, e]) =>
        burden(e.reduce((sum, e) => e.exposure + sum, 0))
      );
      tmbs.sort((a, b) => a - b);
      return { cancer, tmbs };
    })
    .sort((a, b) => (a.cancer > b.cancer ? 1 : b.cancer > a.cancer ? -1 : 0));

  console.log(cancerBurden);
  const totalCancer = cancerBurden.length;
  console.log(totalCancer);

  const absYValue = cancerBurden
    .map((o) => o.tmbs.map((e) => Math.abs(e)))
    .flat();
  const yMax = Math.max(...absYValue);
  const traces = cancerBurden.map((element, index, array) => ({
    element: element,
    index: index,
    array: array,
    name: `${element.cancer}`,
    type: 'scatter',
    marker: { symbol: 'circle-open', size: 3, color: 'black' },
    mode: 'markers',
    y: element.tmbs.map((e) => e),
    average: average(element.tmbs.map((e) => e)),
    hovertemplate: 'Number of mutations: %{y}<br>',
    // x: element.tmbs.map(
    //   (e, i) =>
    //     array
    //       .slice(0, index)
    //       .reduce((x0, curr) => x0 + curr.tmbs.length, 0) +
    //     i +
    //     0.5
    // ),

    x: element.tmbs.map(
      (e, i) => index + 0.1 + (0.8 / element.tmbs.length) * i
    ),
  }));

  console.log('traces:--');
  console.log(traces);

  const topLabel = cancerBurden.map((element, index, array) => ({
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    // x:
    //   array
    //     .slice(0, index)
    //     .reduce((x0, curr) => x0 + curr.tmbs.length, 0) +
    //   (element.tmbs.length - 1) * 0.5,
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

  const bottoLabel1 = cancerBurden.map((element, index, array) => ({
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: (index + index + 1) * 0.5,
    y: -0.1,
    text: `${element.tmbs.length}`,
    showarrow: false,
    font: {
      size: 12,
      color: 'blue',
    },
    align: 'center',
  }));

  const bottoLabelline = cancerBurden.map((element, index, array) => ({
    type: 'line',
    xref: 'x',
    yref: 'paper',
    x0: index + 0.3,
    x1: index + 0.7,
    y0: -0.11,
    y1: -0.11,
    line: {
      width: 1,
      color: 'black',
    },
  }));

  const bottoLabel2 = cancerBurden.map((element, index, array) => ({
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: (index + index + 1) * 0.5,
    y: -0.21,
    text: `${element.tmbs.length}`,
    showarrow: false,
    font: {
      size: 12,
      color: 'green',
    },
    align: 'center',
  }));

  const shapes = cancerBurden.map((element, index, array) => ({
    type: 'rect',
    xref: 'x',
    yref: 'paper',

    // x0: array
    //   .slice(0, index)
    //   .reduce((x0, curr) => x0 + curr.tmbs.length, 0),

    // x1: array
    //   .slice(0, index + 1)
    //   .reduce((x0, curr) => x0 + curr.tmbs.length, 0),
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

  const lines = cancerBurden.map((element, index, array) => ({
    type: 'line',
    xref: 'x',
    yref: 'y',

    // x0: array
    //   .slice(0, index)
    //   .reduce((x0, curr) => x0 + curr.tmbs.length, 0),

    // x1: array
    //   .slice(0, index + 1)
    //   .reduce((x0, curr) => x0 + curr.tmbs.length, 0),
    x0: index + 0.1,
    x1: index + 0.9,

    y0: average(element.tmbs.map((e) => e)),
    y1: average(element.tmbs.map((e) => e)),
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
      autorange: false,
      range: [-Math.ceil(yMax), Math.ceil(yMax)],
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
