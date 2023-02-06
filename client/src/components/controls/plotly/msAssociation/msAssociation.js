import { groupByCustom, linearRegression, round } from '../../utils/utils';
export default function MsAssociation(data, arg) {
  console.log(data);
  console.log(arg);
  const [signatureName1, signatureName2] = arg.signatureName.split(';');
  console.log(signatureName1);
  console.log(signatureName2);
  const groupBySignatureName = groupByCustom(data, (e) => e.signatureName);
  const signatureName1data = groupBySignatureName.get(signatureName1);
  const signatureName2data = groupBySignatureName.get(signatureName2);

  console.log(groupBySignatureName);
  console.log(signatureName1data);
  console.log(signatureName2data);

  function checkBothCal(e, checked) {
    let checkCalculation;
    if (checked) {
      checkCalculation = Math.log10(e['exposure']);
    } else {
      checkCalculation = Math.log10(e['exposure'] + 1);
    }
    return checkCalculation;
  }

  const minX = Math.min(
    ...signatureName1data.map((e) => Math.log10(e['exposure'] + 1))
  );

  const maxX = Math.max(
    ...signatureName1data.map((e) => Math.log10(e['exposure'] + 1))
  );

  const traceSig1 = {
    x: signatureName1data.map((e) => Math.log10(e['exposure'] + 1)),
    name: signatureName1,
    type: 'histogram',
    histnorm: 'density',
    nbinsx: signatureName1data.length - 1,
    yaxis: 'y2',
    marker: { color: '#019E72', line: { color: 'black', width: 1 } },
    hovertemplate:
      '<b>' +
      signatureName1 +
      '</b><br><b>x-Range of ' +
      signatureName1 +
      ' (log10)</b>: %{x}<br><b>Value (log10): </b> %{y}<extra></extra>',
  };

  const traceSig2 = {
    y: signatureName2data.map((e) => Math.log10(e['exposure'] + 1)),
    name: signatureName2,
    type: 'histogram',
    histnorm: 'density',
    nbinsy: signatureName2data.length - 1,
    xaxis: 'x2',
    marker: { color: '#D55E00', line: { color: 'black', width: 1 } },
    hovertemplate:
      '<b>' +
      signatureName2 +
      '</b> <br> <b>x-range of ' +
      signatureName2 +
      ' (log10)</b>: %{y}<br><b>Value (log10): </b> %{x}<extra></extra>',
  };

  const traceMain = {
    x: signatureName1data.map((e) => Math.log10(e['exposure'] + 1)),
    y: signatureName2data.map((e) => Math.log10(e['exposure'] + 1)),
    mode: 'markers',
    type: 'scatter',
    marker: {
      color: '#A3A3A3',
      size: 10,
    },
    opacity: 0.9,
    showlegend: false,
    hovertemplate:
      '<b>Number of mutation in ' +
      signatureName1 +
      ' (log10)</b>: %{x}<br><b>Number of mutation in ' +
      signatureName2 +
      ': (log10)</b> %{y}<extra></extra>',
  };

  const lr = linearRegression(traceMain.x, traceMain.y);
  console.log(lr);

  const traceLine = {
    x: [minX, maxX],
    y: [minX * lr.sl + lr.off, maxX * lr.sl + lr.off],
    name: 'y=' + lr.sl + ' * x + ' + lr.off,
    mode: 'lines',
    marker: {
      color: 'blue',
    },
    hovertemplate:
      '<b>x: </b> %{x}<br><b>y: </b>%{y}<br>' +
      'y=' +
      round(lr.sl, 2) +
      'x + ' +
      round(lr.off, 2) +
      '<extra></extra>',
    showlegend: false,
  };

  const traces = [traceMain, traceLine, traceSig1, traceSig2];

  const detailAnnotation = {
    xref: 'x',
    yref: 'paper',
    x: 0,
    xanchor: 'bottom',
    y: 1,
    yanchor: 'bottom',
    text: 'n<sub>pairs</sub> = ' + signatureName1data.length,
    showarrow: false,
    font: {
      size: 16,
      family: 'Times New Roman',
    },
    align: 'center',
  };

  const layout = {
    showlegend: true,
    hoverlabel: { bgcolor: '#FFF' },
    height: 900,
    bargap: 0,
    autosize: true,
    title: {
      text: '<b>Mutational Signature Association</b>',
    },
    xaxis: {
      anchor: 'y',
      domain: [0.0, 0.83],
      showgrid: true,
      title: {
        text: '<b>Number of mutations in ' + signatureName1 + ' (log10)</b>',
      },
    },
    yaxis: {
      anchor: 'x',
      domain: [0.0, 0.83],
      title: {
        text: '<b>Number of mutations in ' + signatureName2 + ' (log10)</b>',
      },
      showgrid: true,
    },

    xaxis2: { anchor: 'y', domain: [0.85, 1], zerolinecolor: '#EBEBEB' },
    yaxis2: { anchor: 'x', domain: [0.85, 1], zerolinecolor: '#EBEBEB' },

    annotations: [detailAnnotation],
  };
  return { traces: traces, layout: layout };
}
