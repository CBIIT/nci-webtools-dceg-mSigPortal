import { faCommentDollar } from '@fortawesome/free-solid-svg-icons';

export default function TMB(data, tmbTabName, signatureName) {
  function average(arr) {
    const sum = arr.reduce((a, b) => a + b, 0);
    return sum / arr.length || 0;
  }

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
    // average: average(element.samples.map((e) => e.tmb)),
    hovertemplate:
      tmbTabName === 'TMBSignature'
        ? `<b>Sample Name:</b> %{customdata.sampleName}<br><b>Signature Name</b>: ${element.signatureName}` +
          '<br><b>Mutation #: </b>%{y}<extra></extra>'
        : `<b>Sample Name:</b> %{customdata.sampleName}<br><b>Cancer Name:</b>${element.cancer}` +
          '<br><b>Mutation #: </b>%{y}<extra></extra>',
    x: element.samples.map(
      (e, i) => index + 0.1 + (0.8 / element.samples.length) * i
    ),
  }));

  const topLabel = data.map((element, index, array) => ({
    element: element,
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    //x: array.length > 1 ? index : (index + index + 1) * 0.5,
    x: index,
    y: 1.01,
    text:
      tmbTabName === 'TMBSignature' && element.signatureName.length < 13
        ? element.signatureName
        : tmbTabName === 'TMBSignature' && element.signatureName.length > 13
        ? element.signatureName.substring(0, 13) + '...'
        : tmbTabName !== 'TMBSignature' && element.cancer.length < 13
        ? element.cancer
        : element.cancer.substring(0, 13) + '...',
    showarrow: false,
    font: {
      //size: 12,
    },
    align: 'center',
    textangle: 60,
  }));

  const bottoLabel1 = data.map((element, index, array) => ({
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: (index + index + 1) * 0.5,
    y: -0.07,
    text: element.samples.filter(function (x) {
      return x.burden != null;
    }).length,
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
    y0: -0.07,
    y1: -0.07,
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
    y: -0.14,
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
        ...topLabel,
        ...bottoLabel1,
        ...bottoLabel2,
        signatureNameAnnotation,
      ])
    : (annotations = [...topLabel, ...bottoLabel1, ...bottoLabel2]);

  const layout = {
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
    margin: {
      t: 150,
    },
    shapes: [...shapes, ...lines, ...bottoLabelline],
    annotations: annotations,
  };

  //console.log('layout:');
  //console.log(layout);

  var config = {
    //responsive: true,
  };

  return { traces: [...traces], layout: layout, config };
  //return { traces, layout };
}
