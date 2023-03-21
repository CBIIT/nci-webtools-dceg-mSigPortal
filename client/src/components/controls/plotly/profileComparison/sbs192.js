import {
  getTotalMutations,
  getMaxMutations,
  groupDataByMutation,
  getRss,
  getCosineSimilarity,
} from './profileComparison.js';

export default function sbs192(data1, data2) {
  const colors = {
    'C>A': '#03BCEE',
    'C>G': 'black',
    'C>T': '#E32926',
    'T>A': '#CAC9C9',
    'T>C': '#A1CE63',
    'T>G': '#EBC6C4',
  };

  const mutationRegex = /\[(.*)\]/;
  const mutationTypeSort = (a, b) => {
    const mutationTypeRegex = /^\w\:(.*)/;
    return a.mutationType
      .match(mutationTypeRegex)[1]
      .localeCompare(b.mutationType.match(mutationTypeRegex)[1]);
  };
  const formatMutationLabels = (e) => `<b>${e.mutation}</b>`;
  const formatTickLabels = (mutationGroups) =>
    mutationGroups
      .map(({ mutation, data }) =>
        data.map((e) => {
          const color = colors[mutation];
          const regex = /^\w\:(.)\[(.).{2}\](.)$/;
          const match = e.mutationType.match(regex);
          return `${match[1]}<span style="color:${color}"><b>${match[2]}</b></span>${match[3]}`;
        })
      )
      .flat();

  const sample1 = data1[0].sample || data1[0].signatureName;
  const sample2 = data2[0].sample || data2[0].signatureName;

  const mutationGroupSort = (a, b) => {
    const order = Object.keys(colors);
    return order.indexOf(a.mutation) - order.indexOf(b.mutation);
  };

  // get total mutations per sample
  const totalMutations1 = getTotalMutations(data1);
  const totalMutations2 = getTotalMutations(data2);

  // get max mutations per sample
  const maxMutation1 = getMaxMutations(data1) / totalMutations1;
  const maxMutation2 = getMaxMutations(data2) / totalMutations2;
  const maxMutations = Math.max(maxMutation1, maxMutation2);

  // normalize mutations per sample
  const normalizedSample1 = data1.map((e) => ({
    ...e,
    ...(e.mutations >= 0 && {
      mutations: e.mutations / totalMutations1,
    }),
    ...(e.contribution >= 0 && {
      contribution: e.contribution / totalMutations1,
    }),
  }));
  const normalizedSample2 = data2.map((e) => ({
    ...e,
    ...(e.mutations >= 0 && {
      mutations: e.mutations / totalMutations2,
    }),
    ...(e.contribution >= 0 && {
      contribution: e.contribution / totalMutations2,
    }),
  }));

  // separate transcribed and unstranscribed data
  const transcribed1 = normalizedSample1.filter((e) =>
    /^T:/.test(e.mutationType)
  );
  const untranscribed1 = normalizedSample1.filter((e) =>
    /^U:/.test(e.mutationType)
  );
  const transcribed2 = normalizedSample2.filter((e) =>
    /^T:/.test(e.mutationType)
  );
  const untranscribed2 = normalizedSample2.filter((e) =>
    /^U:/.test(e.mutationType)
  );

  const transcribedGroups1 = groupDataByMutation(
    transcribed1,
    mutationRegex,
    mutationGroupSort,
    mutationTypeSort
  );
  const untranscribedGroups1 = groupDataByMutation(
    untranscribed1,
    mutationRegex,
    mutationGroupSort,
    mutationTypeSort
  );
  const transcribedGroups2 = groupDataByMutation(
    transcribed2,
    mutationRegex,
    mutationGroupSort,
    mutationTypeSort
  );
  const untranscribedGroups2 = groupDataByMutation(
    untranscribed2,
    mutationRegex,
    mutationGroupSort,
    mutationTypeSort
  );

  // calcualte difference between samples
  const transcribedDifference = transcribed1.map((e, i) => ({
    ...e,
    ...(e.mutations >= 0 && {
      mutations: e.mutations - transcribed2[i].mutations,
    }),
    ...(e.contribution >= 0 && {
      contribution: e.contribution - transcribed2[i].contribution,
    }),
  }));
  const untranscribedDifference = untranscribed1.map((e, i) => ({
    ...e,
    ...(e.mutations >= 0 && {
      mutations: e.mutations - untranscribed2[i].mutations,
    }),
    ...(e.contribution >= 0 && {
      contribution: e.contribution - untranscribed2[i].contribution,
    }),
  }));

  const transcribedDifferenceGroup = groupDataByMutation(
    transcribedDifference,
    mutationRegex,
    mutationGroupSort,
    mutationTypeSort
  );
  const untranscribedDifferenceGroup = groupDataByMutation(
    untranscribedDifference,
    mutationRegex,
    mutationGroupSort,
    mutationTypeSort
  );

  const transcribedTrace1 = {
    name: 'Transcribed Strand',
    legendgroup: 'transcribed',
    type: 'bar',
    marker: { color: '#004765' },
    x: transcribedGroups1
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
    y: transcribedGroups1
      .map((group, i) => group.data.map((e) => e.mutations || e.contribution))
      .flat(),
    hovertemplate: '<b>Transcribed Strand</b><br> %{x}, %{y} <extra></extra>',
    showlegend: true,
    yaxis: 'y3',
  };
  const untranscribedTrace1 = {
    name: 'Untranscribed Strand',
    legendgroup: 'untranscribed',
    type: 'bar',
    marker: { color: '#E32925' },
    x: untranscribedGroups1
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
    y: untranscribedGroups1
      .map((group, i) => group.data.map((e) => e.mutations || e.contribution))
      .flat(),
    hovertemplate: '<b>Untranscribed Strand</b><br> %{x}, %{y} <extra></extra>',
    showlegend: true,
    yaxis: 'y3',
  };

  const transcribedTrace2 = {
    name: 'Transcribed Strand',
    legendgroup: 'transcribed',
    type: 'bar',
    marker: { color: '#004765' },
    x: transcribedGroups2
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
    y: transcribedGroups2
      .map((group, i) => group.data.map((e) => e.mutations || e.contribution))
      .flat(),
    hovertemplate: '<b>Transcribed Strand</b><br> %{x}, %{y} <extra></extra>',
    showlegend: false,
    yaxis: 'y2',
  };
  const untranscribedTrace2 = {
    name: 'Untranscribed Strand',
    legendgroup: 'untranscribed',
    type: 'bar',
    marker: { color: '#E32925' },
    x: untranscribedGroups2
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
    y: untranscribedGroups2
      .map((group, i) => group.data.map((e) => e.mutations || e.contribution))
      .flat(),
    hovertemplate: '<b>Untranscribed Strand</b><br> %{x}, %{y} <extra></extra>',
    showlegend: false,
    yaxis: 'y2',
  };

  const transcribedDiffernceTrace = {
    name: 'Transcribed Strand',
    legendgroup: 'transcribed',
    type: 'bar',
    marker: { color: '#004765' },
    x: transcribedDifferenceGroup
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
    y: transcribedDifferenceGroup
      .map((group, i) => group.data.map((e) => e.mutations || e.contribution))
      .flat(),
    hovertemplate: '<b>Transcribed Strand</b><br> %{x}, %{y} <extra></extra>',
    showlegend: false,
  };
  const untranscribedDifferenceTrace = {
    name: 'Untranscribed Strand',
    legendgroup: 'untranscribed',
    type: 'bar',
    marker: { color: '#E32925' },
    x: untranscribedDifferenceGroup
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
    y: untranscribedDifferenceGroup
      .map((group, i) => group.data.map((e) => e.mutations || e.contribution))
      .flat(),
    hovertemplate: '<b>Untranscribed Strand</b><br> %{x}, %{y} <extra></extra>',
    showlegend: false,
  };
  const traces = [
    transcribedDiffernceTrace,
    untranscribedDifferenceTrace,
    transcribedTrace2,
    untranscribedTrace2,
    transcribedTrace1,
    untranscribedTrace1,
  ];

  const rss = getRss([...transcribedDifference, ...untranscribedDifference]);
  const cosineSimilarity = getCosineSimilarity(
    normalizedSample1,
    normalizedSample2
  );

  const tickLabels = formatTickLabels(transcribedGroups1);

  const mutationLabelBox = transcribedGroups1.map(
    (group, groupIndex, array) => ({
      type: 'rect',
      xref: 'x',
      yref: 'paper',
      x0: array
        .slice(0, groupIndex)
        .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
      x1: array
        .slice(0, groupIndex + 1)
        .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
      y0: 1.0,
      y1: 1.05,
      fillcolor: colors[group.mutation],
      line: {
        width: 1,
      },
    })
  );
  const sampleBorder1 = transcribedGroups1.map((group, groupIndex, array) => ({
    type: 'rect',
    xref: 'x',
    yref: 'paper',
    x0: array
      .slice(0, groupIndex)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
    x1: array
      .slice(0, groupIndex + 1)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
    y0: 0.67,
    y1: 1,

    line: {
      width: 1,
    },
  }));

  const sampleBorder2 = transcribedGroups2.map((group, groupIndex, array) => ({
    type: 'rect',
    xref: 'x',
    yref: 'paper',
    x0: array
      .slice(0, groupIndex)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
    x1: array
      .slice(0, groupIndex + 1)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
    y0: 0.34,
    y1: 0.66,

    line: {
      width: 1,
    },
  }));

  const differenceBorder = transcribedDifferenceGroup.map(
    (group, groupIndex, array) => ({
      type: 'rect',
      xref: 'x',
      yref: 'paper',
      x0: array
        .slice(0, groupIndex)
        .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
      x1: array
        .slice(0, groupIndex + 1)
        .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
      y0: 0.33,
      y1: 0,

      line: {
        width: 1,
      },
    })
  );
  const differencLabelBox1 = {
    type: 'rect',
    xref: 'paper',
    yref: 'paper',
    x0: 1,
    x1: 1.02,
    y0: 0.67,
    y1: 1,
    fillcolor: '#F0F0F0',
  };
  const sampleLabelBox2 = {
    type: 'rect',
    xref: 'paper',
    yref: 'paper',
    x0: 1,
    x1: 1.02,
    y0: 0.34,
    y1: 0.66,
    fillcolor: '#F0F0F0',
  };

  const differenceLabelBox = {
    type: 'rect',
    xref: 'paper',
    yref: 'paper',
    x0: 1,
    x1: 1.02,
    y0: 0,
    y1: 0.33,
    fillcolor: '#F0F0F0',
  };

  const sampleLabel1 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'middle',
    yanchor: 'middle',
    align: 'center',
    x: 1.0175,
    y: 0.835,
    text: sample1.length > 16 ? sample1.substring(0, 16) + '...' : sample1,
    textangle: 90,
    showarrow: false,
  };

  const sampleLabel2 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'middle',
    yanchor: 'middle',
    align: 'center',
    x: 1.0175,
    y: 0.505,
    text: sample2.length > 16 ? sample2.substring(0, 16) + '...' : sample2,
    textangle: 90,
    showarrow: false,
  };

  const differenceLabel = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'middle',
    yanchor: 'middle',
    align: 'center',
    x: 1.0175,
    y: 0.165,
    text: 'Difference',
    textangle: 90,
    showarrow: false,
  };

  const mutationAnnotation = transcribedGroups1.map(
    (group, groupIndex, array) => ({
      xref: 'x',
      yref: 'paper',
      xanchor: 'bottom',
      yanchor: 'bottom',
      x:
        array
          .slice(0, groupIndex)
          .reduce((lastIndex, b) => lastIndex + b.data.length, 0) +
        (group.data.length - 1) * 0.5,
      y: 1.005,
      text: formatMutationLabels(group),
      showarrow: false,
      font: { color: 'white' },
      align: 'center',
    })
  );

  const yTitleAnnotation = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'middle',
    yanchor: 'middle',
    align: 'center',
    x: -0.051,
    y: 0.5,
    text: '<b>Relative contribution</b>',
    textangle: -90,
    showarrow: false,
  };

  const layout = {
    height: 700,
    hoverlabel: { bgcolor: '#FFF' },
    grid: {
      rows: 3,
      column: 1,
    },
    title:
      '<b>RSS = ' + rss + '; Cosine Similarity = ' + cosineSimilarity + '</b>',
    xaxis: {
      showline: true,
      tickangle: -90,
      tickfont: { family: 'Courier New, monospace' },
      tickmode: 'array',
      tickvals: tickLabels.map((_, i) => i),
      ticktext: tickLabels.map((e) => e),
      linecolor: '#E0E0E0',
      linewidth: 1,
      mirror: 'all',
      ticks: '',
    },
    yaxis: {
      autorange: true,
      linecolor: '#D3D3D3',
      linewidth: 1,
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      showgrid: true,
      gridcolor: '#F5F5F5',
      domain: [0, 0.33],
    },
    yaxis2: {
      autorange: false,
      range: [0, maxMutations * 1.2],
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      domain: [0.34, 0.66],
    },
    yaxis3: {
      autorange: false,
      range: [0, maxMutations * 1.2],
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      domain: [0.67, 1],
    },
    shapes: [
      ...mutationLabelBox,
      differencLabelBox1,
      sampleLabelBox2,
      differenceLabelBox,
      ...sampleBorder1,
      ...sampleBorder2,
      ...differenceBorder,
    ],
    annotations: [
      ...mutationAnnotation,
      sampleLabel1,
      sampleLabel2,
      differenceLabel,
      yTitleAnnotation,
    ],
  };

  return {
    traces,
    layout,
  };
}
