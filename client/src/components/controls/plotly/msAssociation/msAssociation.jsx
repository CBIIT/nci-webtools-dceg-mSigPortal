import {
  linearRegression,
  round,
  calculatePearson,
  calculateSpearman,
} from '@/components/controls/utils/utils';
import groupBy from 'lodash/groupBy';

function formatStatistic(value, decimalPlaces = 2) {
  if (value === Infinity) return '∞';
  if (value === -Infinity) return '−∞';
  return Number.isFinite(value) ? round(value, decimalPlaces) : 'N/A';
}

function formatPValue(value) {
  if (!Number.isFinite(value)) return 'N/A';
  if (value === 0) return '&lt; 0.001';
  return value < 0.001 ? value.toExponential(2) : round(value, 3);
}

function formatCorrelation(label, stats, approximate = false) {
  const pairs = `n<sub>pairs</sub> = ${stats.n}`;
  if (stats.reason) return `${label}: N/A (${stats.reason}), ${pairs}`;
  const notes = [stats.testReason, stats.ciReason].filter(Boolean);
  return (
    `${label}:\tt<sub>Student</sub> = ${formatStatistic(stats.statistic)}` +
    `, p${approximate ? ' (approx.)' : ''} = ${formatPValue(stats.pValue)}` +
    `, r<sub>${label}</sub> = ${formatStatistic(stats.pcorr)}` +
    `, CI<sub>95%${approximate ? ', approx.' : ''}</sub>[` +
    `${formatStatistic(stats.ci[0])}, ${formatStatistic(
      stats.ci[1]
    )}], ${pairs}` +
    (notes.length ? ` (${notes.join('; ')})` : '')
  );
}

export default function MsAssociation(
  data,
  signatureName1,
  signatureName2,
  both = false
) {
  let groupBySample;
  let xValues = [];
  let yValues = [];

  if (both) {
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

  if (xValues.length !== yValues.length) {
    throw new Error(
      `Signature sample counts do not match: "${signatureName1}" has ${xValues.length} samples, "${signatureName2}" has ${yValues.length} samples.`
    );
  }

  const minX = xValues.length
    ? Math.min(...xValues.map((e) => Math.log10(e['exposure'] + 1)))
    : 0;

  const maxX = xValues.length
    ? Math.max(...xValues.map((e) => Math.log10(e['exposure'] + 1)))
    : 0;

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

  const pearsonV = calculatePearson(traceMain.x, traceMain.y);
  const spearman = calculateSpearman(traceMain.x, traceMain.y);
  const lr = pearsonV.reason
    ? null
    : linearRegression(traceMain.x, traceMain.y);
  const traceLine =
    lr && Number.isFinite(lr.sl) && Number.isFinite(lr.off)
      ? {
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
        }
      : null;
  const traces = [traceMain, traceLine, traceSig1, traceSig2].filter(Boolean);

  const detailAnnotation = {
    xref: 'paper',
    yref: 'paper',
    x: 0,
    xanchor: 'bottom',
    y: 1.01,
    yanchor: 'bottom',
    text:
      formatCorrelation('Pearson', pearsonV) +
      '<br>' +
      formatCorrelation(
        'Spearman',
        {
          ...spearman,
          pcorr: spearman.rho,
          statistic: spearman.t,
          ci: [spearman.CILower, spearman.CIUpper],
        },
        true
      ),
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
