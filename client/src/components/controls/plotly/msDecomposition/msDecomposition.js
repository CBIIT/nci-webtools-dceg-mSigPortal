import { groupBy } from 'lodash';

export default function MsDecomposition(data, arg) {
  console.log(data);
  console.log(arg);

  const result = Object.values(data)[0];
  console.log(result);

  const groupByCancer = groupBy(result, 'cancer');
  console.log(groupByCancer);

  function dotp(x, y) {
    function dotp_sum(a, b) {
      return a + b;
    }
    function dotp_times(a, i) {
      return x[i] * y[i];
    }
    return x.map(dotp_times).reduce(dotp_sum, 0);
  }

  function cosineSimilarity(A, B) {
    var similarity =
      dotp(A, B) / (Math.sqrt(dotp(A, A)) * Math.sqrt(dotp(B, B)));
    return similarity;
  }

  function calculate_simularities(
    original_genomes,
    signature,
    signature_activities
  ) {}
  const resultSortedCosine = result.sort((a, b) =>
    a['Cosine_similarity'] > b['Cosine_similarity'] ? 1 : -1
  );
  console.log(resultSortedCosine);

  var trace1 = {
    y: resultSortedCosine.map((e, i) => i),
    x: resultSortedCosine.map((e) => e['Cosine_similarity']),
    fillcolor: '#2B9089',
    line: { color: '#182B2A' },
    fill: 'tozeroy',
    type: 'scatter',
  };

  console.log(trace1);
  const resultSortedL1 = result.sort((a, b) =>
    a['100-L1_Norm_%'] > b['100-L1_Norm_%'] ? 1 : -1
  );
  console.log(resultSortedL1);
  var trace2 = {
    y: resultSortedL1.map((e, i) => i),
    x: resultSortedL1.map((e) => e['100-L1_Norm_%']),
    fill: 'tozeroy',
    type: 'scatter',
    fillcolor: '#2B9089',
    line: { color: '#182B2A' },
    mode: 'lines',
    xaxis: 'x2',
    yaxis: 'y2',
  };
  console.log(trace2);

  const resultSortedL2 = result.sort((a, b) =>
    a['100-L2_Norm_%'] > b['100-L2_Norm_%'] ? 1 : -1
  );
  console.log(resultSortedL2);
  var trace3 = {
    y: resultSortedL2.map((e, i) => i),
    x: resultSortedL2.map((e) => e['100-L2_Norm_%']),
    fill: 'tozeroy',
    type: 'scatter',
    mode: 'lines',
    fillcolor: '#2B9089',
    line: { color: '#182B2A' },
    xaxis: 'x3',
    yaxis: 'y3',
  };

  const resultSortedKL = result.sort((a, b) =>
    a['KL_Divergence'] < b['KL_Divergence'] ? 1 : -1
  );
  console.log(resultSortedKL);
  var trace4 = {
    y: resultSortedKL.map((e, i) => i),
    x: resultSortedKL.map((e) => e['KL_Divergence']),
    fill: 'tonexty',
    type: 'scatter',
    mode: 'lines',
    fillcolor: '#2B9089',
    line: { color: '#182B2A' },
    xaxis: 'x4',
    yaxis: 'y4',
  };
  const resultSortedCorrelation = result.sort((a, b) =>
    a['Correlation'] > b['Correlation'] ? 1 : -1
  );
  console.log(resultSortedKL);
  var trace5 = {
    y: resultSortedCorrelation.map((e, i) => i),
    x: resultSortedCorrelation.map((e) => e['Correlation']),
    fill: 'tonexty',
    type: 'scatter',
    mode: 'lines',
    fillcolor: '#2B9089',
    line: { color: '#182B2A' },
    xaxis: 'x5',
    yaxis: 'y5',
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
  ];

  const boxes = [
    {
      type: 'rect',
      xref: 'paper',
      yref: 'paper',
      x0: 0,
      x1: 0.19,
      y0: 1,
      y1: 1.1,
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
      y1: 1.1,
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
      y1: 1.1,
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
      y1: 1.1,
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
      y1: 1.1,
      fillcolor: '#A1D99B',
      line: {
        width: 0,
      },
    },
  ];

  const traces = [trace1, trace2, trace3, trace4, trace5];

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: 450,
    // width: 1080,
    autosize: true,
    title: {
      text: '<b>Evaluating the Performance of Mutational Signature Decomposition</b>',
    },
    showlegend: false,
    annotations: annotations,
    shapes: boxes,
    xaxis: { domain: [0, 0.19] },
    yaxis2: { anchor: 'x2', showticklabels: false },
    xaxis2: { domain: [0.2, 0.39] },
    yaxis3: { anchor: 'x3', showticklabels: false },
    xaxis3: { domain: [0.4, 0.59] },
    yaxis4: { anchor: 'x4', showticklabels: false },
    xaxis4: { domain: [0.6, 0.79] },
    yaxis5: { anchor: 'x5', showticklabels: false },
    xaxis5: { domain: [0.8, 0.99] },
  };

  return { traces, layout };
}
