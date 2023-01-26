import { groupBy } from 'lodash';
import { groupByCustom } from '../../utils/utils';

export default function MsDecomposition(data, arg) {
  console.log(data);
  console.log(arg);

  const result = Object.values(data)[0];
  console.log(result);

  const grouped = groupByCustom(result, (e) => e.name);
  console.log(grouped);
  console.log(grouped.get('Cosine_similarity'));
  const cosine_similarity = grouped
    .get('Cosine_similarity')
    .sort((a, b) => (a['value'] > b['value'] ? 1 : -1));

  const L1_Norm = grouped
    .get('100-L1_Norm_%')
    .sort((a, b) => (a['value'] > b['value'] ? 1 : -1));

  const L2_Norm = grouped
    .get('100-L2_Norm_%')
    .sort((a, b) => (a['value'] > b['value'] ? 1 : -1));

  const KL_Divergence = grouped
    .get('KL_Divergence')
    .sort((a, b) => (a['value'] < b['value'] ? 1 : -1));

  const Correlation = grouped
    .get('Correlation')
    .sort((a, b) => (a['value'] > b['value'] ? 1 : -1));

  const groupByCancer = groupBy(result, 'cancer');
  console.log(groupByCancer);
  console.log(Object.values(groupByCancer));
  const cancerName =
    Object.keys(groupByCancer) + ' (' + cosine_similarity.length + ')';

  const groupBySample = groupBy(result, 'sample');
  console.log(groupBySample);

  var trace1 = {
    name: 'Cosine Similarity',
    y: cosine_similarity.map((e, i) => i),
    x: cosine_similarity.map((e) => e['value']),
    fillcolor: '#2B9089',
    line: { color: '#182B2A', shape: 'spline', smoothing: 1.3 },
    fill: 'tozeroy',
    type: 'scatter',
    mode: 'lines',
    hovertemplate: '<b>Cosine Similarity:</b> %{x}</b> <extra></extra>',
  };

  console.log(trace1);

  var trace2 = {
    name: '100-L1_Norm_%',
    y: L1_Norm.map((e, i) => i),
    x: L1_Norm.map((e) => e['value']),
    fill: 'tozeroy',
    type: 'scatter',
    line: { color: '#182B2A', shape: 'spline', smoothing: 1.3 },
    mode: 'lines',
    fillcolor: '#2B9089',
    xaxis: 'x2',
    yaxis: 'y',
    hovertemplate: '<b>100-L1_Norm_%:</b> %{x}</b> <extra></extra>',
  };
  console.log(trace2);

  var trace3 = {
    name: '100-L2_Norm_%',
    y: L2_Norm.map((e, i) => i),
    x: L2_Norm.map((e) => e['value']),
    fill: 'tozeroy',
    type: 'scatter',
    line: { color: '#182B2A', shape: 'spline', smoothing: 1.3 },
    mode: 'lines',
    fillcolor: '#2B9089',
    xaxis: 'x3',
    yaxis: 'y',
    hovertemplate: '<b>100-L2_Norm_%: </b>%{x}</b> <extra></extra>',
  };

  var trace4 = {
    name: 'KL_Divergence',
    y: KL_Divergence.map((e, i) => i),
    x: KL_Divergence.map((e) => e['value']),
    fill: 'tozeroy',
    type: 'scatter',
    line: { color: '#182B2A', shape: 'spline', smoothing: 1.3 },
    mode: 'lines',
    fillcolor: '#2B9089',
    xaxis: 'x4',
    yaxis: 'y',
    hovertemplate: '<b>KL_Divergence:</b> %{x}</b> <extra></extra>',
  };

  var trace5 = {
    name: 'Correlation',
    y: Correlation.map((e, i) => i),
    x: Correlation.map((e) => e['value']),
    fill: 'tozeroy',
    type: 'scatter',
    line: { color: '#182B2A', shape: 'spline', smoothing: 1.3 },
    mode: 'lines',
    fillcolor: '#2B9089',
    xaxis: 'x5',
    yaxis: 'y',
    hovertemplate: '<b>Correlation:</b> %{x}</b> <extra></extra>',
  };

  const annotations = [
    {
      xref: 'paper',
      yref: 'paper',
      x: 0.095,
      xanchor: 'center',
      y: 1.02,
      yanchor: 'bottom',
      text: 'Cosine_similarity',
      align: 'center',
      showarrow: false,
    },
    {
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      x: 0.295,
      y: 1.02,
      yanchor: 'bottom',
      text: '100-L1_Norm_%',
      align: 'center',
      showarrow: false,
    },
    {
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      x: 0.495,
      y: 1.02,
      yanchor: 'bottom',
      text: '100-L2_Norm_%',
      align: 'center',
      showarrow: false,
    },
    {
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      x: 0.695,
      y: 1.02,
      yanchor: 'bottom',
      text: 'KL_Divergence',
      align: 'center',
      showarrow: false,
    },
    {
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      x: 0.895,
      y: 1.02,
      yanchor: 'bottom',
      text: 'Correlation',
      align: 'center',
      showarrow: false,
    },

    {
      xref: 'paper',
      yref: 'paper',
      xanchor: 'right',
      x: -0.03,
      y: -0.1,
      yanchor: 'right',
      text: cancerName,
      align: 'center',
      showarrow: false,
    },
  ];

  const boxes = [
    {
      type: 'rect',
      xref: 'paper',
      yref: 'paper',
      x0: 0,
      x1: 0.19,
      y0: 1,
      y1: 1.15,
      fillcolor: '#A1D99B',
      line: {
        width: 0,
      },
    },
    {
      type: 'rect',
      xref: 'paper',
      yref: 'paper',
      x0: 0.2,
      x1: 0.39,
      y0: 1,
      y1: 1.15,
      fillcolor: '#A1D99B',
      line: {
        width: 0,
      },
    },
    {
      type: 'rect',
      xref: 'paper',
      yref: 'paper',
      x0: 0.4,
      x1: 0.59,
      y0: 1,
      y1: 1.15,
      fillcolor: '#A1D99B',
      line: {
        width: 0,
      },
    },
    {
      type: 'rect',
      xref: 'paper',
      yref: 'paper',
      x0: 0.6,
      x1: 0.79,
      y0: 1,
      y1: 1.15,
      fillcolor: '#A1D99B',
      line: {
        width: 0,
      },
    },
    {
      type: 'rect',
      xref: 'paper',
      yref: 'paper',
      x0: 0.8,
      x1: 0.99,
      y0: 1,
      y1: 1.15,
      fillcolor: '#A1D99B',
      line: {
        width: 0,
      },
    },
  ];

  const traces = [trace1, trace2, trace3, trace4, trace5];

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: 350,

    autosize: true,
    title: {
      text: '<b>Evaluating the Performance of Mutational Signature Decomposition</b>',
    },
    showlegend: false,
    annotations: annotations,
    shapes: boxes,
    margin: {
      l: cancerName.length * 10,
    },
    xaxis: { domain: [0, 0.19] },
    yaxis: { anchor: 'x', showticklabels: false, showgrid: false },
    yaxis2: { anchor: 'x2', showticklabels: false, showgrid: false },
    xaxis2: { domain: [0.2, 0.39] },
    yaxis3: { anchor: 'x3', showticklabels: false, showgrid: false },
    xaxis3: { domain: [0.4, 0.59] },
    yaxis4: { anchor: 'x4', showticklabels: false, showgrid: false },
    xaxis4: { domain: [0.6, 0.79] },
    yaxis5: { anchor: 'x5', showticklabels: false, showgrid: false },
    xaxis5: { domain: [0.8, 0.99] },
  };
  console.log(layout);
  return { traces, layout };
}
