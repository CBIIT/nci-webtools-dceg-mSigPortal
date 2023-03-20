import {
  createSampleAnnotation,
  getMaxMutations,
  groupDataByMutation,
  createMutationShapes,
  createMutationAnnotations,
  getTotalMutations,
} from './utils.js';
import { sbsColor } from '../../utils/colors.js';

export default function SBS96(apiData, title = '') {
  const colors = sbsColor;
  const mutationRegex = /\[(.*)\]/;

  const mutationGroupSort = (a, b) => {
    const order = Object.keys(colors);
    return order.indexOf(a.mutation) - order.indexOf(b.mutation);
  };

  const data = groupDataByMutation(apiData, mutationRegex, mutationGroupSort);
  const maxMutation = getMaxMutations(apiData);
  const totalMutations = getTotalMutations(apiData);
  const mutationTypeNames = data
    .map((group) =>
      group.data.map((e) => ({
        mutation: group.mutation,
        mutationType: e.mutationType,
      }))
    )
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
  const mutationAnnotation = createMutationAnnotations(data);
  const mutationShapes = createMutationShapes(data, colors);

  function formatTickLabel(mutation, mutationType) {
    const color = colors[mutation];
    const regex = /^(.)\[(.).{2}\](.)$/;
    const match = mutationType.match(regex);
    return `${match[1]}<span style="color:${color}"><b>${match[2]}</b></span>${match[3]}`;
  }

  const layout = {
    title: `<b>${title}</b>`,
    hoverlabel: { bgcolor: '#FFF' },
    bargap: 0.3,
    height: 450,
    // width: 1080,
    autosize: true,

    xaxis: {
      showline: true,
      tickangle: -90,
      tickfont: {
        family: 'Courier New, monospace',
        color: '#A0A0A0',
      },
      tickmode: 'array',
      tickvals: mutationTypeNames.map((_, i) => i),
      ticktext: mutationTypeNames.map((e) =>
        formatTickLabel(e.mutation, e.mutationType)
      ),
      linecolor: '#E0E0E0',
      linewidth: 1,
      mirror: 'all',
    },
    yaxis: {
      title: {
        text:
          parseFloat(totalMutations).toFixed(2) > 1
            ? '<b>Number of Single Base Substitutions</b>'
            : '<b>Percentage of Single Base Substitutions</b>',
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
      tickformat: parseFloat(totalMutations).toFixed(2) > 1 ? '~s' : '.1%',
      showgrid: true,
      gridcolor: '#F5F5F5',
    },

    shapes: mutationShapes,
    annotations: [...mutationAnnotation, sampleAnnotation],
  };

  return { traces, layout };
}
