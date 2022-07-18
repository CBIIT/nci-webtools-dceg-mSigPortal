export default function DBS78(data, sample) {
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

  const totalMutations = data.reduce(
    (total, mutation) =>
      total +
      mutation.data.reduce((mutationSum, e) => mutationSum + e.mutations, 0),
    0
  );
  const maxMutation = Math.max(
    ...data.map((mutation) => mutation.data.map((e) => e.mutations)).flat()
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
    y: group.data.map((e) => e.mutations),
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

  const sampleAnnotation = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: 0,
    y: 0.88,
    text:
      '<b>' +
      sample +
      ': ' +
      totalMutations.toLocaleString(undefined) +
      ' subs </b>',
    showarrow: false,
    font: {
      size: 18,
    },
    align: 'center',
  };

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: 450,
    //width:1080,
    autosize: true,
    xaxis: {
      showline: true,
      tickangle: -90,
      tickfont: { size: 11 },
      tickmode: 'array',
      tickvals: mutationTypeNames.map((_, i) => i),
      ticktext: mutationTypeNames.map((e) => e.mutationType),
      linecolor: '#E0E0E0',
      linewidth: 1,
      mirror: 'all',
    },
    yaxis: {
      title: 'Number of Double Base Substitutions',
      autorange: false,
      range: [0, maxMutation * 1.2],
      linecolor: '#E0E0E0',
      linewidth: 1,
      tickformat: '~s',
      ticks: 'inside',
      showgrid: true,
      mirror: 'all',
    },

    shapes: shapes,
    annotations: [...mutationAnnotation, sampleAnnotation],
  };

  //console.log("layout");
  //console.log(layout);

  return { traces, layout };
}
