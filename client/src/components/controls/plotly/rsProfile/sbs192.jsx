export default function SBS192(data, sample, tab) {
  const colors = {
    'C>A': '#03BCEE',
    'C>G': 'black',
    'C>T': '#E32926',
    'T>A': '#CAC9C9',
    'T>C': '#A1CE63',
    'T>G': '#EBC6C4',
  };

  const transcribed = data.filter((e) => /^T:/.test(e.mutationType));
  const untranscribed = data.filter((e) => /^U:/.test(e.mutationType));

  const totalMutations =
    transcribed.reduce((total, e) => total + e.mutations, 0) +
    untranscribed.reduce((total, e) => total + e.mutations, 0);

  const maxMutation = Math.max(
    ...[
      ...transcribed.map((e) => e.mutations),
      ...untranscribed.map((e) => e.mutations),
    ]
  );

  const maxVal = Math.max(...data.map((o) => o?.mutations || o?.contribution));

  const groupByMutationT = transcribed.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.mutationType.match(mutationRegex)[1];
    const signature = {
      mutationType: e.mutationType.substring(2, e.mutationType.length),
      contribution: e?.mutations || e?.contribution,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];

    return groups;
  }, {});
  const groupByMutationU = untranscribed.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.mutationType.match(mutationRegex)[1];
    const signature = {
      mutationType: e.mutationType.substring(2, e.mutationType.length),
      contribution: e?.mutations || e?.contribution,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  const transformU = Object.entries(groupByMutationU).map(
    ([mutation, data]) => ({
      mutation,
      data,
    })
  );

  const flatSortedT = Object.values(groupByMutationT)
    .flat()
    .sort((a, b) =>
      a.mutationType
        .substring(0, 5)
        .localeCompare(b.mutationType.substring(0, 5))
    )
    .sort((a, b) =>
      a.mutationType
        .substring(2, a.mutationType.length)
        .localeCompare(b.mutationType.substring(2, a.mutationType.length))
    )
    .sort((a, b) =>
      a.mutationType
        .substring(0, 5)
        .localeCompare(b.mutationType.substring(0, 5))
    )
    .sort((a, b) =>
      a.mutationType
        .substring(2, 5)
        .localeCompare(b.mutationType.substring(2, 5))
    );
  const flatSortedU = Object.values(groupByMutationU)
    .flat()
    .sort((a, b) =>
      a.mutationType
        .substring(0, 5)
        .localeCompare(b.mutationType.substring(0, 5))
    )
    .sort((a, b) =>
      a.mutationType
        .substring(2, 5)
        .localeCompare(b.mutationType.substring(2, 5))
    );

  const tracesT = {
    name: 'Transcrribed Strand',
    type: 'bar',
    marker: { color: '#004765' },
    x: flatSortedT.map((element, index, array) => index),
    y: flatSortedT.map((element, index, array) => element.contribution),
    hovertemplate: '<b>Transcribed Strand</b><br> %{x}, %{y}<extra></extra>',
    //hoverinfo: 'x+y',
    showlegend: true,
  };
  const tracesU = {
    name: 'Untranscribed Strand',
    type: 'bar',
    marker: { color: '#E32925' },
    x: flatSortedU.map((element, index, array) => index),
    y: flatSortedU.map((element, index, array) => element.contribution),
    hovertemplate: '<b>Transcribed Strand</b><br> %{x}, %{y} <extra></extra>',
    //hoverinfo: 'x+y',
    showlegend: true,
  };
  const annotations = Object.entries(groupByMutationU).map(
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
      text: `<b>${mutation}</b>`,
      showarrow: false,
      font: {
        size: 18,
      },
      align: 'center',
    })
  );

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
          ' transcribed subs</b>',
    showarrow: false,
    font: {
      size: 24,
      family: 'Arial',
    },
    align: 'center',
  };

  const shapes1 = Object.entries(groupByMutationU).map(
    ([mutation, _], groupIndex, array) => ({
      type: 'rect',
      xref: 'x',
      yref: 'paper',
      x0: array
        .slice(0, groupIndex)
        .reduce((x0, [_, sigs]) => x0 + sigs.length, -0.4),
      x1: array
        .slice(0, groupIndex + 1)
        .reduce((x0, [_, sigs]) => x0 + sigs.length, -0.6),
      y0: 1.05,
      y1: 1.01,
      fillcolor: colors[mutation],
      line: {
        width: 0,
      },
    })
  );

  const shapes2 = Object.entries(groupByMutationU).map(
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
      fillcolor: colors[mutation],
      line: {
        width: 0,
      },
      opacity: 0.15,
    })
  );
  const traces = [tracesT, tracesU];

  function formatTickLabel(mutation, mutationType) {
    const color = colors[mutation];
    const regex = /^(.)\[(.).{2}\](.)$/;
    const match = mutationType.match(regex);
    return `${match[1]}<span style="color:${color}"><b>${match[2]}</b></span>${match[3]}`;
  }

  const mutationTypeNames = transformU
    .map((group) =>
      group.data.map((e) => ({
        mutation: group.mutation,
        mutationType: e.mutationType,
      }))
    )
    .flat();

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    showlegend: true,
    height: 450,
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
        color: '#A0A0A0',
      },
      tickmode: 'array',
      // tickvals: [...flatSortedU.map((_, i) => i)],
      // ticktext: [...flatSortedU.map((e) => e.mutationType)],
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
        text: '<b>Percent of Single Base Substitutions</b>',
        font: {
          family: 'Times New Roman',
        },
      },
      autorange: false,
      range: [0, maxVal + maxVal * 0.2],
      linecolor: '#E0E0E0',
      linewidth: 1,
      mirror: 'all',
      tickformat: '.1%',
    },

    shapes: [...shapes1, ...shapes2],
    annotations: [...annotations, sampleAnnotation],
  };

  return { traces, layout };
}
