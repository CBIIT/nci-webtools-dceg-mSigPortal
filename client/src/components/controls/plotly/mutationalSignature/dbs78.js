export default function DBS78(data, sample) {
  const colors = {
    'AC>': '#09BCED',
    'AT>': '#0266CA',
    'CC>': '#9FCE62',
    'CG>': '#006501',
    'CT>': '#FF9898',
    'GC>': '#E22925',
    'TA>': '#FEB065',
    'TC>': '#FD8000',
    'TG>': '#CB98FD',
    'TT>': '#4C0299',
  };

  const totalMutations = data.reduce((a, e) => a + parseInt(e.mutations), 0);
  const maxVal = Math.max(...data.map((o) => o.mutations));
  //console.log(maxVal);
  const numberWithCommas = (x) =>
    x.toString().replace(/\B(?<!\.\d*)(?=(\d{3})+(?!\d))/g, ',');
  //console.log(totalMutations);
  //console.log("data");
  //console.log(data);

  const dataFilter = data.filter(
    (item) =>
      item.mutationType.substring(0, 2) === 'AC' ||
      item.mutationType.substring(0, 2) === 'AT' ||
      item.mutationType.substring(0, 2) === 'CC' ||
      item.mutationType.substring(0, 2) === 'CG' ||
      item.mutationType.substring(0, 2) === 'CT' ||
      item.mutationType.substring(0, 2) === 'GC' ||
      item.mutationType.substring(0, 2) === 'TA' ||
      item.mutationType.substring(0, 2) === 'TC' ||
      item.mutationType.substring(0, 2) === 'TG' ||
      item.mutationType.substring(0, 2) === 'TT'
  );

  //console.log(dataFilter);
  // group data by dominant mutation
  const groupByMutation = dataFilter.reduce((groups, e, i) => {
    const mutation = e.mutationType.substring(0, 3);
    const signature = {
      mutationType: e.mutationType,
      contribution: e.mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  const flatSorted = Object.values(groupByMutation).flat();

  const traces = Object.entries(groupByMutation).map(
    ([mutation, signatures], groupIndex, array) => ({
      name: mutation,
      type: 'bar',
      marker: { color: colors[mutation] },
      //x: signatures.map((e) => e.mutationType),
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

  //console.log(traces);
  const annotations = Object.entries(groupByMutation).map(
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
      text: `<b>${mutation}NN</b>`,
      showarrow: false,
      font: {
        size: 18,
      },
      align: 'center',
    })
  );
  //console.log(annotations);
  const shapes = Object.entries(groupByMutation).map(
    ([mutation, signatures], groupIndex, array) => ({
      type: 'rect',
      xref: 'x',
      yref: 'paper',
      x0: array
        .slice(0, groupIndex)
        .reduce((x0, [_, signatures]) => x0 + signatures.length, -0.4),
      // x0: groupIndex * 16 - 0.4,
      x1: array
        .slice(0, groupIndex + 1)
        .reduce((x0, [_, signatures]) => x0 + signatures.length, -0.6),
      // x1: groupIndex * 16 + signatures.length - 0.6,
      y0: 1.05,
      y1: 1.01,
      fillcolor: colors[mutation],
      line: {
        width: 0,
      },
    })
  );
  //console.log(shapes);

  const sampleAnnotation = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: 0,
    y: 0.92,
    text:
      '<b>' +
      sample +
      ': ' +
      numberWithCommas(totalMutations) +
      ' double subs</b>',
    showarrow: false,
    font: {
      size: 18,
    },
    align: 'center',
  };

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: 500,
    width: 1080,
    xaxis: {
      //title: "Double Substitution",
      showticklabels: true,
      showline: true,
      tickangle: -90,
      tickfont: {
        size: 10,
      },
      tickmode: 'array',
      tickvals: flatSorted.map((_, i) => i),
      ticktext: flatSorted.map((e) =>
        e.mutationType.substring(
          e.mutationType.length - 2,
          e.mutationType.length
        )
      ),
      linecolor: 'black',
      linewidth: 1,
      mirror: true,
    },
    yaxis: {
      title: 'Number of Double Base Substitutions',
      autorange: false,
      range: [0, maxVal + maxVal * 0.2],
      linecolor: 'black',
      linewidth: 1,
      mirror: true,
    },

    shapes: shapes,
    annotations: [...annotations, sampleAnnotation],
  };

  //console.log("layout");
  //console.log(layout);

  return { traces, layout };
}
