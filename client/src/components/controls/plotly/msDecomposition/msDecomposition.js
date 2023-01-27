import { groupBy } from 'lodash';
import { groupByCustom } from '../../utils/utils';

export default function MsDecomposition(data, arg) {
  console.log(data);
  console.log(arg);

  const result = Object.values(data)[0];
  console.log(result);

  const grouped = groupByCustom(result, (e) => e.name);
  console.log(grouped);

  const cosine_similarity = grouped.get('Cosine_similarity');
  //.sort((a, b) => (a['value'] > b['value'] ? 1 : -1));
  const L1_Norm = grouped.get('100-L1_Norm_%');
  //.sort((a, b) => (a['value'] > b['value'] ? 1 : -1));
  const L2_Norm = grouped.get('100-L2_Norm_%');
  //  .sort((a, b) => (a['value'] > b['value'] ? 1 : -1));

  const KL_Divergence = grouped.get('KL_Divergence');
  // .sort((a, b) => (a['value'] > b['value'] ? 1 : -1));

  const Correlation = grouped.get('Correlation');
  // .sort((a, b) => (a['value'] > b['value'] ? 1 : -1));

  const groupByCancer = groupBy(result, 'cancer');

  const cancerName =
    Object.keys(groupByCancer) + ' (' + Object.values(data)[1].length + ')';

  const traces = [
    {
      type: 'violin',
      spanmode: 'soft',
      hoveron: 'points+kde',
      // span: [0, 5],
      side: 'positive', //positive side means right for vertical violin plots
      //y: 4,
      x: cosine_similarity.map((e) => e['value']),
      points: 'none',
      box: {
        visible: false,
      },
      boxpoints: false,
      line: {
        color: 'black',
      },
      fillcolor: '#2B9089',
      scalemode: 'count',
      marker: {
        line: {
          width: 0,
        },
        symbol: 'line-ns',
      },
      opacity: 0.6,
      meanline: {
        visible: true,
      },
      y0: 'Cosine_Similarity',
    },
    {
      type: 'violin',
      spanmode: 'soft',
      hoveron: 'points+kde',
      //span: [0,5],
      side: 'positive', //positive side means right for vertical violin plots
      //y: 4,
      x: L1_Norm.map((e) => e['value']),
      points: 'none',
      box: {
        visible: false,
      },
      boxpoints: false,
      line: {
        color: 'black',
      },
      fillcolor: '#8491B4',
      scalemode: 'count',
      marker: {
        line: {
          width: 0,
        },
        symbol: 'line-ns',
      },
      opacity: 0.6,
      meanline: {
        visible: true,
      },
      y0: '100-L1_Norm%',
      xaxis: 'x2',
      yaxis: 'y2',
    },
    {
      type: 'violin',
      spanmode: 'soft',
      hoveron: 'points+kde',
      //span: [0,5],
      side: 'positive', //positive side means right for vertical violin plots
      //y: 4,
      x: L2_Norm.map((e) => e['value']),
      points: 'none',
      box: {
        visible: false,
      },
      boxpoints: false,
      line: {
        color: 'black',
      },
      fillcolor: '#B24745',
      scalemode: 'count',
      marker: {
        line: {
          width: 0,
        },
        symbol: 'line-ns',
      },
      opacity: 0.6,
      meanline: {
        visible: true,
      },
      y0: '100-L2_Norm%',
      xaxis: 'x3',
      yaxis: 'y3',
    },
    {
      type: 'violin',
      spanmode: 'soft',
      hoveron: 'points+kde',
      //span: [0,5],
      side: 'positive', //positive side means right for vertical violin plots
      //y: 4,
      x: KL_Divergence.map((e) => e['value']),
      points: 'none',
      box: {
        visible: false,
      },
      boxpoints: false,
      line: {
        color: 'black',
      },
      fillcolor: '#E4AF69',
      scalemode: 'count',
      marker: {
        line: {
          width: 0,
        },
        symbol: 'line-ns',
      },
      opacity: 0.6,
      meanline: {
        visible: true,
      },
      y0: 'KL_Divergence',
      xaxis: 'x4',
      yaxis: 'y4',
    },
    {
      type: 'violin',
      spanmode: 'soft',
      hoveron: 'points+kde',
      //span: [0,5],
      side: 'positive', //positive side means right for vertical violin plots
      //y: 4,
      x: Correlation.map((e) => e['value']),
      points: 'none',
      box: {
        visible: false,
      },
      boxpoints: false,
      line: {
        color: 'black',
      },
      fillcolor: '#AE1F63',
      scalemode: 'count',
      marker: {
        line: {
          width: 0,
        },
        symbol: 'line-ns',
      },
      opacity: 0.6,
      meanline: {
        visible: true,
      },
      y0: 'Correlation',
      xaxis: 'x5',
      yaxis: 'y5',
    },
  ];

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
    yaxis: { anchor: 'x', showticklabels: false },
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
