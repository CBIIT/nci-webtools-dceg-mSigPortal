export default function SBS96(data, sample) {
  const colors = {
    'C>A': '#03BCEE',
    'C>G': 'black',
    'C>T': '#E32926',
    'T>A': '#CAC9C9',
    'T>C': '#A1CE63',
    'T>G': '#EBC6C4',
  };

  const numberWithCommas = (x) =>
    x.toString().replace(/\B(?<!\.\d*)(?=(\d{3})+(?!\d))/g, ',');

  const totalMutations = data.reduce(
    (total, base) =>
      total +
      base.mutationTypes.reduce((baseSum, e) => baseSum + e.mutations, 0),
    0
  );
  const maxMutation = Math.max(
    ...data.map((base) => base.mutationTypes.map((e) => e.mutations)).flat()
  );

  const mutationTypeNames = data
    .map((group) =>
      group.mutationTypes.map((e) => ({
        base: group.base,
        mutationType: e.mutationType,
      }))
    )
    .flat();

  const traces = data.map((group, groupIndex, array) => ({
    name: group.base,
    type: 'bar',
    marker: { color: colors[group.base] },
    x: [...group.mutationTypes.keys()].map(
      (e) =>
        e +
        array
          .slice(0, groupIndex)
          .reduce((lastIndex, b) => lastIndex + b.mutationTypes.length, 0)
    ),
    y: group.mutationTypes.map((e) => e.mutations),
    hoverinfo: 'x+y',
    showlegend: false,
  }));

  const baseAnnotation = data.map((group, groupIndex, array) => ({
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x:
      array
        .slice(0, groupIndex)
        .reduce((lastIndex, b) => lastIndex + b.mutationTypes.length, 0) +
      (group.mutationTypes.length - 1) * 0.5,
    y: 1.04,
    text: `<b>${group.base}</b>`,
    showarrow: false,
    font: { size: 18 },
    align: 'center',
  }));

  const sampleAnnotation = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: 0,
    y: 0.88,
    text:
      '<b>' + sample + ': ' + numberWithCommas(totalMutations) + ' subs </b>',
    showarrow: false,
    font: {
      size: 18,
    },
    align: 'center',
  };

  const shapes = data.map((group, groupIndex, array) => ({
    type: 'rect',
    xref: 'x',
    yref: 'paper',
    x0: array
      .slice(0, groupIndex)
      .reduce((lastIndex, e) => lastIndex + e.mutationTypes.length, -0.35),
    x1: array
      .slice(0, groupIndex + 1)
      .reduce((lastIndex, e) => lastIndex + e.mutationTypes.length, -0.65),
    y0: 1.05,
    y1: 1.01,
    fillcolor: colors[group.base],
    line: {
      width: 0,
    },
  }));

  function formatTickLabel(base, mutationType) {
    const color = colors[base];
    const regex = /^(.)\[(.).{2}\](.)$/;
    const match = mutationType.match(regex);
    return `${match[1]}<span style="color:${color}">${match[2]}</span>${match[3]}`;
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
      tickfont: { size: 11 },
      tickmode: 'array',
      tickvals: mutationTypeNames.map((_, i) => i),
      ticktext: mutationTypeNames.map((e) =>
        formatTickLabel(e.base, e.mutationType)
      ),
      linecolor: '#E0E0E0',
      linewidth: 1,
      mirror: 'all',
    },
    yaxis: {
      title: 'Number of Single Base Substitutions',
      autorange: false,
      range: [0, maxMutation + maxMutation * 0.2],
      ticks: 'inside',
      linecolor: '#E0E0E0',
      linewidth: 1,
      mirror: 'all',
      tickformat: '~s',
      showgrid: false,
    },

    shapes: shapes,
    annotations: [...baseAnnotation, sampleAnnotation],
  };

  return { traces, layout };
}
