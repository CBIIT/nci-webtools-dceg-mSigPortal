import { groupBy } from 'lodash';

export default function MsDecomposition(data, arg) {
  console.log(data);
  console.log(arg);

  const result = Object.values(data)[0];
  console.log(result);

  const groupByCancer = groupBy(result, 'cancer');
  console.log(groupByCancer);
  console.log(Object.values(groupByCancer));
  const cancerName =
    Object.keys(groupByCancer) +
    ' (' +
    Object.values(groupByCancer)[0].length +
    ')';
  console.log(cancerName);
  console.log(cancerName.length);
  const resultSortedCosine = result.sort((a, b) =>
    a['Cosine_similarity'] > b['Cosine_similarity'] ? 1 : -1
  );
  console.log(resultSortedCosine);

  var trace1 = {
    name: 'Cosine Similarity',
    y: resultSortedCosine.map((e, i) => i),
    x: resultSortedCosine.map((e) => e['Cosine_similarity']),
    fillcolor: '#2B9089',
    line: { color: '#182B2A', shape: 'spline', smoothing: 1.3 },
    fill: 'tozeroy',
    type: 'scatter',
    mode: 'lines',
    hovertemplate: '<b>Cosine Similarity:</b> %{x}</b> <extra></extra>',
  };

  console.log(trace1);
  const resultSortedL1 = result.sort((a, b) =>
    a['100-L1_Norm_%'] > b['100-L1_Norm_%'] ? 1 : -1
  );
  console.log(resultSortedL1);
  var trace2 = {
    name: '00-L1_Norm_%',
    y: resultSortedL1.map((e, i) => i),
    x: resultSortedL1.map((e) => e['100-L1_Norm_%']),
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

  const resultSortedL2 = result.sort((a, b) =>
    a['100-L2_Norm_%'] > b['100-L2_Norm_%'] ? 1 : -1
  );
  console.log(resultSortedL2);
  var trace3 = {
    name: '100-L2_Norm_%',
    y: resultSortedL2.map((e, i) => i),
    x: resultSortedL2.map((e) => e['100-L2_Norm_%']),
    fill: 'tozeroy',
    type: 'scatter',
    line: { color: '#182B2A', shape: 'spline', smoothing: 1.3 },
    mode: 'lines',
    fillcolor: '#2B9089',
    xaxis: 'x3',
    yaxis: 'y',
    hovertemplate: '<b>100-L2_Norm_%: </b>%{x}</b> <extra></extra>',
  };

  const resultSortedKL = result.sort((a, b) =>
    a['KL_Divergence'] < b['KL_Divergence'] ? 1 : -1
  );
  console.log(resultSortedKL);
  var trace4 = {
    name: 'KL_Divergence',
    y: resultSortedKL.map((e, i) => i),
    x: resultSortedKL.map((e) => e['KL_Divergence']),
    fill: 'tozeroy',
    type: 'scatter',
    line: { color: '#182B2A', shape: 'spline', smoothing: 1.3 },
    mode: 'lines',
    fillcolor: '#2B9089',
    xaxis: 'x4',
    yaxis: 'y',
    hovertemplate: '<b>KL_Divergence:</b> %{x}</b> <extra></extra>',
  };
  const resultSortedCorrelation = result.sort((a, b) =>
    a['Correlation'] > b['Correlation'] ? 1 : -1
  );
  console.log(resultSortedKL);
  var trace5 = {
    name: 'Correlation',
    y: resultSortedCorrelation.map((e, i) => i),
    x: resultSortedCorrelation.map((e) => e['Correlation']),
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
