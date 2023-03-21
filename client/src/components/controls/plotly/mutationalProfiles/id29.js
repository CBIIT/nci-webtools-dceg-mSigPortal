import {
  createSampleAnnotation,
  getMaxMutations,
  groupDataByMutation,
  createMutationShapes,
  createMutationAnnotations,
} from './utils.js';

export default function ID29(apiData, title = '') {
  const colors = {
    '[+C]': '#0072b2',
    '[+T]': '#d55e00',
    '[+>1]': '#cc79a7',
    '[-C]': '#56b4e9',
    '[-T]': '#e69f00',
    '[->1]': '#009e73',
    '[-]': '#911eb4',
  };
  const mutationTypeSort = (a, b) => {
    const mutationTypeOrder = [
      'A',
      'G',
      'C',
      'T',
      'CC',
      'TT',
      'LR',
      'NonR',
      'Rep',
      'MH',
    ];
    const mutationRegex = /\[.*\](.*)/;

    return (
      mutationTypeOrder.indexOf(a.mutationType.match(mutationRegex)[1]) -
      mutationTypeOrder.indexOf(b.mutationType.match(mutationRegex)[1])
    );
  };
  const groupRegex = /(\[.*\])/;

  const mutationGroupSort = (a, b) => {
    const order = Object.keys(colors);
    return order.indexOf(a.mutation) - order.indexOf(b.mutation);
  };

  const data = groupDataByMutation(
    apiData,
    groupRegex,
    mutationGroupSort,
    mutationTypeSort
  );
  const maxMutation = getMaxMutations(apiData);

  const mutationTypeNames = data
    .map((group) =>
      group.data.map((e) => ({
        mutation: group.mutation,
        mutationType: e.mutationType,
      }))
    )
    .flat();

  const traces = data.map((group, i, array) => ({
    name: group.mutation,
    type: 'bar',
    marker: { color: colors[group.mutation] },
    x: [...group.data.keys()].map(
      (e) =>
        e +
        array.slice(0, i).reduce((lastIndex, b) => lastIndex + b.data.length, 0)
    ),
    y: group.data.map((e) => e.mutations || e.contribution),
    hoverinfo: 'x+y',
    showlegend: false,
  }));

  const sampleAnnotation = createSampleAnnotation(apiData);
  const mutationAnnotation = createMutationAnnotations(data);
  const mutationShapes = createMutationShapes(data, colors);

  const layout = {
    title: `<b>${title}</b>`,
    hoverlabel: { bgcolor: '#FFF' },
    bargap: 0.3,
    height: 450,
    // width: 1080,
    autosize: true,

    xaxis: {
      showline: true,
      tickangle: 45,
      tickfont: {
        family: 'Courier New, monospace',
        // color: '#A0A0A0',
      },
      tickmode: 'array',
      tickvals: mutationTypeNames.map((_, i) => i),
      ticktext: mutationTypeNames.map((e) => `<b>${e.mutationType}</b>`),
      linecolor: '#E0E0E0',
      linewidth: 1,
      mirror: 'all',
    },
    yaxis: {
      title: {
        text: '<b>Mutation Probability</b>',
        font: {
          family: 'Times New Roman',
        },
      },
      autorange: false,
      range: [0, maxMutation * 1.2],
      ticks: 'inside',
      tickcolor: '#D3D3D3',
      linecolor: '#D3D3D3',
      linewidth: 1,
      mirror: 'all',
      tickformat: Number.isInteger(traces[0].y[0]) ? '~s' : '.1%',
      showgrid: true,
      gridcolor: '#F5F5F5',
    },

    shapes: mutationShapes,
    annotations: [...mutationAnnotation, sampleAnnotation],
  };

  return { traces, layout };
}
