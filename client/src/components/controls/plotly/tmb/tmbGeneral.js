export default function TMB(data, tmbTabName, signatureName) {
  const totalCancer = data.length;

  const absYValue = data
    .map((o) => o.samples.map((e) => Math.abs(e.burden)))
    .flat();
  const yMax = Math.max(...absYValue);

  const traces = data.map((element, index, array) => ({
    element: element,
    cancer: element.cancer,
    type: 'scatter',
    marker: { symbol: 'circle-open', size: 4, color: 'black' },
    mode: 'markers',
    customdata: element.samples.map((e) => ({
      sampleName: e.sample,
    })),
    y: element.samples.map((e) => e.burden),
    hovertemplate:
      tmbTabName === 'TMBSignature'
        ? `<b>Sample Name:</b> %{customdata.sampleName}<br><b>Signature Name</b>: ${element.signatureName}` +
          '<br><b>Mutation #: </b>%{y}<extra></extra>'
        : `<b>Sample Name:</b> %{customdata.sampleName}<br><b>Cancer Name:</b>${element.cancer}` +
          '<br><b>Mutation #: </b>%{y}<extra></extra>',
    x: element.samples.map(
      (e, i) => index + 0.1 + (0.8 / element.samples.length) * i
    ),
    showlegend: false,
  }));

  const traceLabelGreen = {
    type: 'bar',
    y: [0],
    x: [0],
    marker: {
      size: 0,
      color: ['green'],
    },
    name: 'Total number of samples',
    hoverinfo: 'skip',
  };

  const groupLabel = data.map((element, index, array) => ({
    element: element,
    xref: 'x',
    yref: 'paper',
    xanchor: 'left',
    yanchor: 'top',
    x: index + 0.2,
    y: 0,
    text:
      tmbTabName === 'TMBSignature' ? element.signatureName : element.cancer,
    showarrow: false,
    font: {
      //size: 12,
    },
    align: 'left',
    textangle: totalCancer > 4 ? 60 : 0,
  }));

  const totalCountAnnotation = data.map((element, index, array) => ({
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: (index + index + 1) * 0.5,
    //y: -0.14,
    y: 1.01,
    text: `${element.totalSamples}`,
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

  const lines = data.map((element, index, array) => ({
    type: 'line',
    xref: 'x',
    yref: 'y',
    x0: index + 0.1,
    x1: index + 0.9,
    y0: element.medianBurden,
    y1: element.medianBurden,
    line: {
      width: 1,
      color: 'red',
    },
  }));

  const signatureNameAnnotation = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: 0.01,
    y: 0.9,
    text: signatureName,
    showarrow: false,
    font: {
      size: 18,
      family: 'Arial',
    },
    align: 'center',
  };
  let annotations = [];
  signatureName != null
    ? (annotations = [
        ...groupLabel,
        ...totalCountAnnotation,
        signatureNameAnnotation,
      ])
    : (annotations = [...groupLabel, ...totalCountAnnotation]);

  // find the longest label to calculate extra height margin
  const labels = groupLabel.map((e) => e.text);
  const longest = labels.reduce((a, e) => (a > e.length ? a : e.length), 0);
  const extraMargin = longest < 10 ? 60 : longest * 7.5;

  const layout = {
    width: totalCancer > 4 ? null : 400 + 150 * (totalCancer - 1),
    autosize: true,
    height: 500 + extraMargin,
    legend: { orientation: 'h', x: totalCancer > 4 ? 0.4 : 0.2, y: 1.15 },
    //title: { text: 'Tumor Mutational Burden' },
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
    margin: {
      b: extraMargin,
      t: 100,
      l: 0,
    },
    shapes: [...shapes, ...lines],
    annotations: annotations,
  };

  var config = {
    //responsive: true,
  };

  return {
    traces: [...traces, traceLabelGreen],
    layout: layout,
    config,
  };
}
