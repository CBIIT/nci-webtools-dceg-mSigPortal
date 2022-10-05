import {
  createSampleAnnotation,
  getMaxMutations,
  groupDataByMutation,
  createMutationShapes,
  createMutationAnnotations,
} from './utils';

export default function DBS78(apiData) {
  const colors = {
    AC: '#09BCED',
    AT: '#0266CA',
    CC: '#9FCE62',
    CG: '#006501',
    CT: '#FF9898',
    GC: '#E22925',
    TA: '#FEB065',
    TC: '#FD8000',
    TG: '#CB98FD',
    TT: '#4C0299',
  };
  const mutationRegex = /^(.{2})/;
  const mutationOrder = Object.keys(colors);

  const data = groupDataByMutation(apiData, mutationRegex, mutationOrder);
  const maxMutation = getMaxMutations(apiData);

  const mutationTypeNames = data
    .map((group) => group.data.map((e) => e.mutationType.slice(-2)))
    .flat();

  const traces = data.map((group, groupIndex, array) => ({
    name: group.mutation,
    type: 'bar',
    marker: { color: colors[group.mutation] },
    x: [...group.data.keys()].map(
      (e) =>
        e +
        array
          .slice(0, groupIndex)
          .reduce((lastIndex, b) => lastIndex + b.data.length, 0)
    ),
    y: group.data.map((e) => e.mutations || e.contribution),
    hoverinfo: 'x+y',
    showlegend: false,
  }));

  const sampleAnnotation = createSampleAnnotation(apiData);
  const mutationAnnotation = createMutationAnnotations(data, '>NN');
  const mutationShapes = createMutationShapes(data, colors);

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: 450,
    //width:1080,
    autosize: true,
    xaxis: {
      showline: true,
      tickangle: -90,
      tickfont: {
        family: 'Courier New, monospace',
        size: 14,
      },
      tickmode: 'array',
      tickvals: mutationTypeNames.map((_, i) => i),
      ticktext: mutationTypeNames.map((e) => e),
      linecolor: 'black',
      linewidth: 1,
      mirror: 'all',
    },
    yaxis: {
      title: {
        text: apiData[0].mutations
          ? '<b>Number of Double Base Substitutions</b>'
          : '<b>Percentage of Double Base Substitutions</b>',
        font: {
          family: 'Times New Roman',
        },
      },
      autorange: false,
      range: [0, maxMutation * 1.2],
      linecolor: 'black',
      linewidth: 1,
      tickformat: apiData[0].contribution ? '.1%' : '~s',
      ticks: 'inside',
      tickcolor: '#D3D3D3',
      showgrid: true,
      mirror: 'all',
      gridcolor: '#F5F5F5',
    },

    shapes: mutationShapes,
    annotations: [...mutationAnnotation, sampleAnnotation],
  };

  return { traces, layout };
}
