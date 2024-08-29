export default function DBS78(rawData, sample) {
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
  const dbsdata = ['AC', 'AT', 'CC', 'CG', 'CT', 'GC', 'TA', 'TC', 'TG', 'TT'];
  const groupByMutation = rawData.reduce((acc, e, i) => {
    const mutationRegex = /^(.{2})/;
    const mutation = e.mutationType.match(mutationRegex)[1];

    acc[mutation] = acc[mutation] ? [...acc[mutation], e] : [e];
    return acc;
  }, {});

  const data = Object.entries(groupByMutation).map(([mutation, data]) => ({
    mutation,
    data,
  }));
  const datafilter = data.filter((e) => dbsdata.includes(e.mutation));
  const totalMutations = datafilter.reduce(
    (total, mutation) =>
      total +
      mutation.data.reduce((mutationSum, e) => mutationSum + e.contribution, 0),
    0
  );
  const maxMutation = Math.max(
    ...datafilter
      .map((mutation) => mutation.data.map((e) => e.contribution))
      .flat()
  );
  const mutationTypeNames = datafilter
    .map((group) =>
      group.data.map((e) => ({
        mutation: group.mutation,
        mutationType: e.mutationType,
      }))
    )
    .flat();

  const traces = datafilter.map((group, groupIndex, array) => ({
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

  const mutationAnnotation = datafilter.map((group, groupIndex, array) => ({
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
    text: `<b>${group.mutation}>NN</b>`,
    showarrow: false,
    font: { size: 18 },
    align: 'center',
  }));

  const shapes = datafilter.map((group, groupIndex, array) => ({
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
    x: 0.01,
    y: 0.9,
    text: '<b>' + sample + ' </b>',
    showarrow: false,
    font: {
      size: 24,
      family: 'Arial',
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
      tickfont: {
        family: 'Courier New, monospace',
        size: 14,
      },
      tickmode: 'array',
      tickvals: mutationTypeNames.map((_, i) => i),
      ticktext: mutationTypeNames.map((e) => e.mutationType.slice(-2)),
      linecolor: 'black',
      linewidth: 1,
      mirror: 'all',
    },
    yaxis: {
      title: {
        text: '<b>Percent of Double Base Substitutions</b>',
        font: {
          family: 'Times New Roman',
        },
      },
      autorange: false,
      range: [0, maxMutation * 1.2],
      linecolor: 'black',
      linewidth: 1,
      //tickformat: maxMutation > 1000 ? '~s' : '',
      tickformat: '.1%',
      ticks: 'inside',
      tickcolor: '#D3D3D3',
      showgrid: true,
      mirror: 'all',
      gridcolor: '#F5F5F5',
    },

    shapes: shapes,
    annotations: [...mutationAnnotation, sampleAnnotation],
  };

  return { traces, layout };
}
