export default function SBS96(rawData, sample, tab) {
  const colors = {
    'C>A': '#03BCEE',
    'C>G': 'black',
    'C>T': '#E32926',
    'T>A': '#CAC9C9',
    'T>C': '#A1CE63',
    'T>G': '#EBC6C4',
  };
  //const { profile, matrix, sample } = args;

  const groupByMutation = rawData.reduce((acc, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.mutationType.match(mutationRegex)[1];

    acc[mutation] = acc[mutation] ? [...acc[mutation], e] : [e];
    return acc;
  }, {});

  const data = Object.entries(groupByMutation).map(([mutation, data]) => ({
    mutation,
    data,
  }));

  const totalMutations = data.reduce(
    (total, mutation) =>
      total +
      mutation.data.reduce((mutationSum, e) => mutationSum + e.contribution, 0),
    0
  );
  const maxMutation = Math.max(
    ...data.map((mutation) => mutation.data.map((e) => e.contribution)).flat()
  );

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
    y: group.data.map((e) => e.contribution),
    hoverinfo: 'x+y',
    showlegend: false,
  }));

  const mutationAnnotation = data.map((group, groupIndex, array) => ({
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x:
      array
        .slice(0, groupIndex)
        .reduce((lastIndex, b) => lastIndex + b.data.length, 0) +
      (group.data.length - 1) * 0.5,
    y: 1.04,
    text: `<b>${group.mutation}</b>`,
    showarrow: false,
    font: { size: 18 },
    align: 'center',
  }));

  const sampleAnnotation = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: 0.01,
    y: 0.88,
    text:
      tab === 'rsProfile'
        ? '<b>' + sample + '</b>'
        : '<b>' +
          sample +
          ': ' +
          totalMutations.toLocaleString(undefined) +
          ' subs </b>',
    showarrow: false,
    font: {
      size: 24,
      family: 'Arial',
    },
    align: 'center',
  };

  const shapes = data.map((group, groupIndex, array) => ({
    type: 'rect',
    xref: 'x',
    yref: 'paper',
    x0: array
      .slice(0, groupIndex)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.35),
    x1: array
      .slice(0, groupIndex + 1)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.65),
    y0: 1.05,
    y1: 1.01,
    fillcolor: colors[group.mutation],
    line: {
      width: 0,
    },
  }));

  function formatTickLabel(mutation, mutationType) {
    const color = colors[mutation];
    const regex = /^(.)\[(.).{2}\](.)$/;
    const match = mutationType.match(regex);
    return `${match[1]}<span style="color:${color}"><b>${match[2]}</b></span>${match[3]}`;
  }

  const layout = {
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
        text: '<b>Percentage of Single Base Substitutions</b>',
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
      tickformat: '.1%',
      showgrid: true,
      gridcolor: '#F5F5F5',
    },

    shapes: shapes,
    annotations: [...mutationAnnotation, sampleAnnotation],
  };

  return { traces, layout };
}
