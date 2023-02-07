export default function SBS288(data, sample, tab) {
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
  const neutral = data.filter((e) => /^N:/.test(e.mutationType));

  const totalMutations =
    transcribed.reduce((total, e) => total + e.contribution, 0) +
    untranscribed.reduce((total, e) => total + e.contribution, 0) +
    neutral.reduce((total, e) => total + e.contribution, 0);
  const maxMutation = Math.max(
    ...[
      ...transcribed.map((e) => e?.mutations || e?.contribution),
      ...untranscribed.map((e) => e?.mutations || e?.contribution),
    ]
  );
  //// ------ bar char left  --------- //

  const groupByMutationWoFirstLetter = data.reduce((groups, e, i) => {
    const mutation = e.mutationType.substring(2, e.mutationType.length);
    const signature = {
      mutationType: e.mutationType,
      contribution: e.contribution,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  const totalMutationsGroup = Object.entries(groupByMutationWoFirstLetter).map(
    ([mutation, signatures], groupIndex, array) => ({
      mutationType: mutation,
      signatures: signatures,
      total: signatures.reduce(
        (a, e) => a + parseInt(e?.mutations || e?.contribution),
        0
      ),
    })
  );

  const groupByTotal = totalMutationsGroup.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.mutationType.match(mutationRegex)[1];
    const signature = {
      mutationType: e.mutationType,
      contribution: e.total,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  const flatSortedTotal = Object.values(groupByTotal).flat();
  const maxValTotal = Math.max(...flatSortedTotal.map((o) => o.contribution));

  const tracesBarTotal = Object.entries(groupByTotal).map(
    ([mutation, signatures], groupIndex, array) => ({
      name: mutation,
      type: 'bar',
      marker: { color: colors[mutation] },
      //   x: signatures.map((e) => e.mutationType),
      //x: signatures.map((e, i) => groupIndex * signatures.length + i),
      x: signatures.map(
        (e, i) =>
          array
            .slice(0, groupIndex)
            .reduce((x0, [_, sigs]) => x0 + sigs.length, 0) + i
      ),
      y: signatures.map((e) => e.contribution),
      hoverinfo: 'x+y',
      showlegend: false,
    })
  );

  //////------------- bar chart right ---------------//////

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

  const totalMutationsGroupT = Object.entries(groupByMutationT).map(
    ([mutation, signatures], groupIndex, array) => ({
      mutationType: mutation,
      signatures: signatures,
      total: signatures.reduce((a, e) => a + parseInt(e.contribution), 0),
    })
  );

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

  const totalMutationsGroupU = Object.entries(groupByMutationU).map(
    ([mutation, signatures], groupIndex, array) => ({
      mutationType: mutation,
      signatures: signatures,
      total: signatures.reduce((a, e) => a + parseInt(e.contribution), 0),
    })
  );

  const groupByMutationN = neutral.reduce((groups, e, i) => {
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

  const totalMutationsGroupN = Object.entries(groupByMutationN).map(
    ([mutation, signatures], groupIndex, array) => ({
      mutationType: mutation,
      signatures: signatures,
      total: signatures.reduce(
        (a, e) => a + parseInt(e?.mutations || e?.contribution),
        0
      ),
    })
  );

  const groupByFirstLetter = data.reduce((groups, e, i) => {
    const mutation = e.mutationType.substring(0, 1);
    const signature = {
      mutationType: e.mutationType,
      contribution: e?.mutations || e?.contribution,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];

    return groups;
  }, {});

  const totalGroupByFirstLetter = Object.entries(groupByFirstLetter).map(
    ([mutation, signatures], groupIndex, array) => ({
      mutationType: mutation,
      signatures: signatures,
      total: signatures.reduce(
        (a, e) => a + parseInt(e?.mutations || e?.contribution),
        0
      ),
    })
  );

  const flatSortedT = Object.values(totalMutationsGroupT).flat();
  const flatSortedU = Object.values(totalMutationsGroupU).flat();
  const flatSortedN = Object.values(totalMutationsGroupN).flat();
  const flatSortedFirstLetter = Object.values(totalGroupByFirstLetter).flat();

  const maxValByTotal = Math.max(
    ...[
      ...flatSortedT.map((e) => e.total),
      ...flatSortedU.map((e) => e.total),
      ...flatSortedN.map((e) => e.total),
      ...flatSortedFirstLetter.map((e) => e.total),
    ]
  );
  Object.entries(flatSortedFirstLetter).forEach(
    ([key, value], groupIndex, array) => {
      if (value.mutationType === 'T') {
        value.mutationType = 'All';
        flatSortedT.push(value);
      } else if (value.mutationType === 'U') {
        value.mutationType = 'All';
        flatSortedU.push(value);
      } else {
        value.mutationType = 'All';
        flatSortedN.push(value);
      }
    }
  );

  const transcribedTraces = {
    name: 'Transcrribed',
    type: 'bar',
    marker: { color: '#004765' },

    x: flatSortedT.map((element, index, array) => element.total),
    y: flatSortedT.map(
      (element, index, array) => `<b>${element.mutationType}<b>`
    ),
    xaxis: 'x2',
    yaxis: 'y2',
    //hoverinfo: 'x2+y2',
    hovertemplate: '<b>Transcrribed</b><br>%{y}, %{x} <extra></extra>',
    showlegend: true,
    orientation: 'h',
  };

  const untranscribedTraces = {
    name: 'Untranscribed',
    type: 'bar',
    marker: { color: '#E32925' },
    x: flatSortedU.map((element, index, array) => element.total),
    y: flatSortedU.map(
      (element, index, array) => `<b>${element.mutationType}<b>`
    ),
    xaxis: 'x2',
    yaxis: 'y2',
    //hoverinfo: 'x2+y2',
    hovertemplate: '<b>Untranscribed</b><br>%{y}, %{x} <extra></extra>',
    showlegend: true,
    orientation: 'h',
  };

  const neutralTraces = {
    name: 'Nontranscribed',
    type: 'bar',
    marker: { color: '#008001' },
    x: flatSortedN.map((element, index, array) => element.total),
    y: flatSortedN.map(
      (element, index, array) => `<b>${element.mutationType}<b>`
    ),
    xaxis: 'x2',
    yaxis: 'y2',
    //hoverinfo: 'x2+y2',
    hovertemplate: '<b>Nontranscribed</b><br>%{y}, %{x} <extra></extra>',
    showlegend: true,
    orientation: 'h',
  };

  const traces = [
    ...tracesBarTotal,
    neutralTraces,
    untranscribedTraces,
    transcribedTraces,
  ];

  const annotations = Object.entries(groupByTotal).map(
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
          ' subs </b>',
    showarrow: false,
    font: {
      size: 24,
      family: 'Arial',
    },
    align: 'center',
  };

  const transformU = Object.entries(groupByMutationU).map(
    ([mutation, data]) => ({
      mutation,
      data,
    })
  );

  const mutationTypeNames = transformU
    .map((group) =>
      group.data.map((e) => ({
        mutation: group.mutation,
        mutationType: e.mutationType,
      }))
    )
    .flat();

  function formatTickLabel(mutation, mutationType) {
    const color = colors[mutation];
    const regex = /^(.)\[(.).{2}\](.)$/;
    const match = mutationType.match(regex);
    return `${match[1]}<span style="color:${color}"><b>${match[2]}</b></span>${match[3]}`;
  }

  const shapes = Object.entries(groupByTotal).map(
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
      mutation: mutation,
    })
  );

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    bargap: 0.3,
    height: 450,
    //width:1080,
    autosize: true,
    legend: {
      x: 1,
      xanchor: 'right',
      y: 0,
      traceorder: 'reversed',
      bgcolor: '#FFFFFF',
      bordercolor: '#E1E1E1',
      borderwidth: 1,
    },
    grid: {
      rows: 1,
      columns: 2,
      pattern: 'independent',
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
      domain: [0, 0.75],
    },
    yaxis: {
      title: {
        text: '<b>Percent of Single Base Substitutions</b>',
        font: {
          family: 'Times New Roman',
        },
      },
      autorange: false,
      range: [0, maxValTotal + maxValTotal * 0.2],
      tickcolor: '#D3D3D3',
      linecolor: '#D3D3D3',
      linewidth: 1,
      tickformat: maxValTotal >= 1000 ? '~s' : '',
      ticks: 'inside',
      showgrid: true,
      gridcolor: '#F5F5F5',
      mirror: 'all',
    },

    xaxis2: {
      showticklabels: true,
      showline: true,
      tickfont: {
        size: 12,
      },
      domain: [0.8, 1],
      //tickformat: maxValByTotal >= 1000 ? '~s' : '',
      tickformat: '.1%',
      ticks: 'outside',
      linecolor: '#E0E0E0',
      linewidth: 1,
      showgrid: false,
    },
    yaxis2: {
      showline: true,
      tickangle: 0,
      tickfont: {
        size: 12,
      },
      anchor: 'x2',
      categoryorder: 'category descending',
      tickformat: '~s',
      ticks: 'outside',
      linecolor: '#D3D3D3',
      tickcolor: '#D3D3D3',
      linewidth: 1,
      showgrid: false,
    },

    shapes: shapes,
    annotations: [...annotations, sampleAnnotation],
  };

  return { traces, layout };
}
