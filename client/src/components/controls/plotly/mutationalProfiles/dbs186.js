import { createSampleAnnotation } from './utils';

export default function DBS186(data, title = '') {
  const colors = {
    'CC>': '#09BCEE',
    'CT>': '#A0CE63',
    'TC>': '#FE9898',
    'TT>': '#FE8002',
  };

  const arrayDataT = [];
  const arrayDataU = [];

  Object.values(data).forEach((group) => {
    if (group.mutationType.substring(0, 1) === 'T') {
      arrayDataT.push(group);
    } else if (group.mutationType.substring(0, 1) === 'U') {
      arrayDataU.push(group);
    }
  });

  const totalMutations = [...arrayDataT, ...arrayDataU].reduce(
    (a, e) => a + parseInt(e.mutations),
    0
  );
  const dataUT = [...arrayDataT, ...arrayDataU];

  const maxVal = Math.max(...dataUT.map((o) => o.mutations));

  // group data by dominant mutation
  const T_groupByMutation = arrayDataT.reduce((groups, e, i) => {
    const mutation = e.mutationType.substring(2, 4);
    const signature = {
      mutationType: e.mutationType,
      contribution: e.mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  const T_flatSorted = Object.values(T_groupByMutation).flat();

  const tracesT = {
    name: 'Transcrribed Strand',
    type: 'bar',
    marker: { color: '#004765' },
    x: T_flatSorted.map((element, index, array) => index),
    y: T_flatSorted.map((element, index, array) => element.contribution),
    hovertemplate: '<b>Transcribed Strand</b><br>%{x}, %{y}<extra></extra>',
    //hoverinfo: 'x+y',
    showlegend: true,
  };

  const U_groupByMutation = arrayDataU.reduce((groups, e, i) => {
    const mutation = e.mutationType.substring(2, 4);
    const signature = {
      mutationType: e.mutationType,
      contribution: e.mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  const U_flatSorted = Object.values(U_groupByMutation).flat();

  const tracesU = {
    name: 'Untranscribed Strand',
    type: 'bar',
    marker: { color: '#E32925' },
    x: U_flatSorted.map((element, index, array) => index),
    y: U_flatSorted.map((element, index, array) => element.contribution),
    hovertemplate: '<b>Untranscribed Strand</b><br>%{x}, %{y}<extra></extra>',
    //hoverinfo: 'x+y',
    showlegend: true,
  };

  const traces = [tracesT, tracesU];

  const annotations = Object.entries(T_groupByMutation).map(
    ([mutation, signatures], groupIndex, array) => ({
      xref: 'x',
      yref: 'paper',
      xanchor: 'bottom',
      yanchor: 'bottom',
      x:
        array
          .slice(0, groupIndex)
          .reduce((x0, [_, sigs]) => x0 + sigs.length, 0) +
        (signatures.length - 1) * 0.5,
      y: 1.04,
      text: `<b>${mutation}>NN</b>`,
      showarrow: false,
      font: {
        size: 18,
      },
      align: 'center',
    })
  );

  const shapes1 = Object.entries(T_groupByMutation).map(
    ([mutation, signatures], groupIndex, array) => ({
      type: 'rect',
      xref: 'x',
      yref: 'paper',
      x0: array
        .slice(0, groupIndex)
        .reduce((x0, [_, signatures]) => x0 + signatures.length, -0.4),
      // x0: groupIndex * 16 - 0.4,
      x1: array
        .slice(0, groupIndex + 1)
        .reduce((x0, [_, signatures]) => x0 + signatures.length, -0.6),
      // x1: groupIndex * 16 + signatures.length - 0.6,
      y0: 1.05,
      y1: 1.01,
      fillcolor:
        colors[
          signatures[0].mutationType.substring(
            2,
            signatures[0].mutationType.length - 2
          )
        ],
      line: {
        width: 0,
      },
    })
  );
  const shapes2 = Object.entries(T_groupByMutation).map(
    ([mutation, signatures], groupIndex, array) => ({
      type: 'rect',
      xref: 'x',
      yref: 'paper',
      y0: 1,
      x0: array
        .slice(0, groupIndex)
        .reduce((x0, [_, sigs]) => x0 + sigs.length, -0.4),
      x1: array
        .slice(0, groupIndex + 1)
        .reduce((x0, [_, sigs]) => x0 + sigs.length, -0.6),
      y1: 0,
      fillcolor:
        colors[
          signatures[0].mutationType.substring(
            2,
            signatures[0].mutationType.length - 2
          )
        ],
      line: {
        width: 0,
      },
      opacity: 0.15,
    })
  );

  const sampleAnnotation = createSampleAnnotation(data);

  const layout = {
    title: `<b>${title}</b>`,
    hoverlabel: { bgcolor: '#FFF' },
    showlegend: true,
    height: 600,
    //width:1080,
    autosize: true,
    legend: {
      x: 1,
      xanchor: 'right',
      y: 1,
      bgcolor: '#FFFFFF',
      bordercolor: '#E1E1E1',
      borderwidth: 1,
    },
    xaxis: {
      showticklabels: true,
      showline: true,
      tickangle: -90,
      tickfont: {
        family: 'Courier New, monospace',
        size: 22,
      },
      tickmode: 'array',
      tickvals: T_flatSorted.map((_, i) => i),
      ticktext: T_flatSorted.map((e) =>
        e.mutationType.substring(
          e.mutationType.length - 2,
          e.mutationType.length
        )
      ),
      linecolor: '#E0E0E0',
      linewidth: 1,
      mirror: 'all',
    },
    yaxis: {
      title: {
        text: '<b>Number of Double Base Substitutions</b>',
        font: {
          family: 'Times New Roman',
          size: 22,
        },
      },
      autorange: false,
      range: [0, maxVal + maxVal * 5],
      ticks: 'inside',
      linecolor: '#E0E0E0',
      tickcolor: '#D3D3D3',
      linewidth: 1,
      mirror: 'all',
      tickfont: {
        size: 16,
      },
      showgrid: true,
      gridcolor: '#F5F5F5',
    },

    shapes: [...shapes1, ...shapes2],
    annotations: [...annotations, sampleAnnotation],
  };

  return { traces, layout };
}
