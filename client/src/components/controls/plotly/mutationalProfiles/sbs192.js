import {
  createSampleAnnotation,
  getMaxMutations,
  groupDataByMutation,
  createMutationShapes,
  createMutationAnnotations,
} from './utils';
import { sbsColor } from '../../utils/colors';

export default function SBS192(apiData, title = '') {
  const colors = sbsColor;
  const mutationRegex = /\[(.*)\]/;
  const mutationTypeSort = (a, b) => {
    const mutationTypeRegex = /^\w\:(.*)/;
    return a.mutationType
      .match(mutationTypeRegex)[1]
      .localeCompare(b.mutationType.match(mutationTypeRegex)[1]);
  };
  const mutationGroupSort = (a, b) => {
    const order = Object.keys(colors);
    return order.indexOf(a.mutation) - order.indexOf(b.mutation);
  };

  const transcribed = apiData.filter((e) => /^T:/.test(e.mutationType));
  const untranscribed = apiData.filter((e) => /^U:/.test(e.mutationType));

  const transcribedGroups = groupDataByMutation(
    transcribed,
    mutationRegex,
    mutationGroupSort,
    mutationTypeSort
  );
  const untranscribedGroups = groupDataByMutation(
    untranscribed,
    mutationRegex,
    mutationGroupSort,
    mutationTypeSort
  );
  const maxMutation = getMaxMutations(apiData);

  const transcribedTraces = {
    name: 'Transcribed Strand',
    type: 'bar',
    marker: { color: '#004765' },
    x: transcribedGroups
      .map((group, i, array) =>
        [...group.data.keys()].map(
          (e) =>
            e +
            array
              .slice(0, i)
              .reduce((lastIndex, b) => lastIndex + b.data.length, 0)
        )
      )
      .flat(),
    y: transcribedGroups
      .map((group, i) => group.data.map((e) => e.mutations || e.contribution))
      .flat(),
    hovertemplate: '<b>Transcribed Strand</b><br> %{x}, %{y} <extra></extra>',
    showlegend: true,
  };
  const untranscribedTraces = {
    name: 'Untranscribed Strand',
    type: 'bar',
    marker: { color: '#E32925' },
    x: untranscribedGroups
      .map((e, i, array) =>
        [...e.data.keys()].map(
          (e) =>
            e +
            array
              .slice(0, i)
              .reduce((lastIndex, b) => lastIndex + b.data.length, 0)
        )
      )
      .flat(),
    y: untranscribedGroups
      .map((group, i) => group.data.map((e) => e.mutations || e.contribution))
      .flat(),
    hovertemplate: '<b>Untranscribed Strand</b><br> %{x}, %{y} <extra></extra>',
    showlegend: true,
  };

  const sampleAnnotation = createSampleAnnotation(
    apiData,
    'Transcribed Substitutions'
  );
  const mutationAnnotation = createMutationAnnotations(transcribedGroups);
  const mutationShapes = createMutationShapes(transcribedGroups, colors);
  const backgroundShapes = mutationShapes.map((e) => ({
    ...e,
    y0: 1,
    y1: 0,
    opacity: 0.15,
  }));

  const traces = [transcribedTraces, untranscribedTraces];

  function formatTickLabel(mutation, mutationType) {
    const color = colors[mutation];
    const regex = /^\w\:(.)\[(.).{2}\](.)$/;
    const match = mutationType.match(regex);
    return `${match[1]}<span style="color:${color}"><b>${match[2]}</b></span>${match[3]}`;
  }

  const mutationTypeNames = transcribedGroups
    .map((group) =>
      group.data.map((e) => ({
        mutation: group.mutation,
        mutationType: e.mutationType,
      }))
    )
    .flat();

  const layout = {
    title: `<b>${title}</b>`,
    hoverlabel: { bgcolor: '#FFF' },
    showlegend: true,
    height: 450,
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
        text: apiData[0].mutations
          ? '<b>Number of Double Base Substitutions</b>'
          : '<b>Percentage of Double Base Substitutions</b>',
        font: {
          family: 'Times New Roman',
        },
      },
      autorange: false,
      range: [0, maxMutation * 1.2],
      linecolor: '#E0E0E0',
      linewidth: 1,
      mirror: 'all',
      tickformat: Number.isInteger(traces[0].y[0]) ? '~s' : '.1%',
    },

    shapes: [...mutationShapes, ...backgroundShapes],
    annotations: [...mutationAnnotation, sampleAnnotation],
  };

  return { traces, layout };
}
