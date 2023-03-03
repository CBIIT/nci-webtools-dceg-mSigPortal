import {
  linearRegression,
  round,
  calculatePearson,
  calculateSpearman,
} from '../../utils/utils';
import { groupBy } from 'lodash';
import pcorrtest from '@stdlib/stats-pcorrtest';

export default function MsAssociation(data, arg) {
  const [signatureName1, signatureName2] = arg.signatureName.split(';');
  const checked = arg.both;

  let groupBySample;
  let xValues = [];
  let yValues = [];

  if (checked) {
    const dataFilter = groupBy(
      data.filter((o) => o['exposure'] > 0),
      'sample'
    );
    if (signatureName1 === signatureName2) {
      groupBySample = { ...dataFilter };
    } else {
      groupBySample = Object.values(dataFilter).filter((e) => e.length > 1);
    }
  } else {
    groupBySample = groupBy(data, 'sample');
  }

  const dataArraySample = Object.values(groupBySample);

  for (var i = 0; i < dataArraySample.length; i++) {
    for (var j = 0; j < dataArraySample[i].length; j++) {
      if (dataArraySample[i][j].signatureName === signatureName1) {
        xValues.push(dataArraySample[i][j]);
      }
      if (dataArraySample[i][j].signatureName === signatureName2) {
        yValues.push(dataArraySample[i][j]);
      }
    }
  }

  const minX = Math.min(...xValues.map((e) => Math.log10(e['exposure'] + 1)));

  const maxX = Math.max(...xValues.map((e) => Math.log10(e['exposure'] + 1)));

  const traceSig1 = {
    //x: signatureName1data.map((e) => Math.log10(e['exposure'] + 1)),
    x: xValues.map((e) => Math.log10(e['exposure'] + 1)),
    name: signatureName1,
    type: 'histogram',
    histnorm: 'density',
    // nbinsx: round(data.length / 1.75),
    nbinsx: 35,
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
    //y: signatureName2data.map((e) => Math.log10(e['exposure'] + 1)),
    y: yValues.map((e) => Math.log10(e['exposure'] + 1)),
    name: signatureName2,
    type: 'histogram',
    histnorm: 'density',
    // nbinsy: round(data.length / 1.75),
    nbinsy: 35,
    xaxis: 'x2',
    marker: { color: '#D55E00', line: { color: 'black', width: 1 } },
    hovertemplate:
      '<b>' +
      signatureName2 +
      '</b> <br> <b>x-range of ' +
      signatureName2 +
      ' (log10)</b> %{y}<br><b>Value (log10): </b> %{x}<extra></extra>',
  };

  const traceMain = {
    x: xValues.map((e) => Math.log10(e['exposure'] + 1)),
    y: yValues.map((e) => Math.log10(e['exposure'] + 1)),
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

  let pearsonV;
  if (traceMain.x.length > 3) {
    pearsonV = pcorrtest(traceMain.x, traceMain.y);
  } else {
    pearsonV = calculatePearson(traceMain.x, traceMain.y);
  }
  const spearman = calculateSpearman(traceMain.x, traceMain.y);
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
    xref: 'paper',
    yref: 'paper',
    x: 0,
    xanchor: 'bottom',
    y: 1.01,
    yanchor: 'bottom',
    text:
      'Pearson:\tt<sub>Student</sub> = ' +
      round(pearsonV.statistic, 2) +
      ', p = ' +
      round(pearsonV.pValue, 3) +
      ', r<sub>Pearson</sub> = ' +
      round(pearsonV.pcorr, 2) +
      ', CI<sub>95%</sub>[' +
      round(pearsonV.ci[0], 2) +
      ', ' +
      round(pearsonV.ci[1], 2) +
      '], n<sub>pairs</sub> = ' +
      dataArraySample.length +
      '<br>Spearman:\tt<sub>Student</sub> =' +
      round(spearman.t, 2) +
      ', p = ' +
      round(spearman.pValue, 3) +
      ', r<sub>Spearman</sub> = ' +
      round(spearman.rho, 2) +
      ', CI<sub>95%</sub>[' +
      round(spearman.CILower, 2) +
      ', ' +
      round(spearman.CIUpper, 2) +
      '], n<sub>pairs</sub> = ' +
      spearman.n,
    showarrow: false,
    font: {
      size: 16,
      family: 'Times New Roman',
    },
    align: 'left',
  };

  const layout = {
    showlegend: true,
    hoverlabel: { bgcolor: '#FFF' },
    height: 700,
    bargap: 0,
    autosize: true,
    title: {
      text: '<b>Mutational Signature Association</b>',
      font: {
        family: 'Arial',
        size: 18,
      },
    },
    legend: {
      title: { text: '\t Signature Names:' },
    },
    xaxis: {
      domain: [0.0, 0.83],

      showgrid: true,
      title: {
        text: '<b>Number of mutations in ' + signatureName1 + ' (log10)</b>',
      },
    },
    yaxis: {
      domain: [0.0, 0.83],
      title: {
        text: '<b>Number of mutations in ' + signatureName2 + ' (log10)</b>',
      },
      showgrid: true,
    },

    xaxis2: { anchor: 'y', domain: [0.85, 1], zerolinecolor: '#EBEBEB' },
    yaxis2: { anchor: 'x', domain: [0.85, 1], zerolinecolor: '#EBEBEB' },

    annotations: [detailAnnotation],
    margin: {
      t: 150,
    },
  };
  return { traces: traces, layout: layout };
}
