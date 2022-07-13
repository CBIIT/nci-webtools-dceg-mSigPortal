export default function SBS96(data, sample) {
  const colors = {
    'C>A': '#03BCEE',
    'C>G': 'black',
    'C>T': '#E32926',
    'T>A': '#CAC9C9',
    'T>C': '#A1CE63',
    'T>G': '#EBC6C4',
  };
  //console.log("data--:");
  //console.log(data);
  const numberWithCommas = (x) =>
    x.toString().replace(/\B(?<!\.\d*)(?=(\d{3})+(?!\d))/g, ',');

  const totalMutations = data.reduce((a, e) => a + parseInt(e.mutations), 0);
  const maxVal = Math.max(...data.map((o) => o.mutations));
  // console.log("maxVal--:");
  // console.log(maxVal);

  // group data by dominant mutation
  const groupByMutation = data.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.mutationType.match(mutationRegex)[1];
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

  // console.log("groupByMutation:");
  // console.log(groupByMutation);
  //console.log("FlatSorted");
  //console.log(flatSorted);
  const traces = Object.entries(groupByMutation).map(
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
  // console.log("traces:");
  // console.log(traces);

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
      text: `<b>${mutation}</b>`,
      showarrow: false,
      font: {
        size: 18,
      },
      align: 'center',
    })
  );

  const xannotations = flatSorted.map((num, index) => ({
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: index,
    y: -0.1,
    text: num.mutationType.replace(
      /\[(.*)\]/,
      num.mutationType.substring(2, 3)
    ),
    showarrow: false,
    font: {
      size: 10,
      color: colors[num.mutationType.substring(2, 5)],
    },
    align: 'center',
    num: num,
    index: index,
    textangle: -90,
  }));

  const sampleAnnotation = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: 0,
    y: 0.92,
    text:
      '<b>' + sample + ': ' + numberWithCommas(totalMutations) + ' subs </b>',
    showarrow: false,
    font: {
      size: 18,
    },
    align: 'center',
  };

  const shapes = Object.entries(groupByMutation).map(
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
  // console.log("shapes");
  // console.log(shapes);

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    bargap: 0.3,
    height: 450,
    width: 1080,

    xaxis: {
      showticklabels: false,
      showline: true,
      tickangle: -90,
      tickfont: {
        size: 10,
      },
      tickmode: 'array',
      tickvals: flatSorted.map((_, i) => i),
      // ticktext: flatSorted.map((e) =>
      //   e.mutationType.replace(/\[(.*)\]/, e.mutationType.substring(2, 3))
      // ),
      ticktext: flatSorted.map((e) => e.mutationType),
      linecolor: 'black',
      linewidth: 1,
      mirror: true,
      tickformat: '~s',
    },
    yaxis: {
      title: 'Number of Single Base Substitutions',
      autorange: false,
      range: [0, maxVal + maxVal * 0.15],
      linecolor: 'black',
      linewidth: 1,
      mirror: true,
      tickformat: '~s',
    },

    shapes: shapes,
    annotations: [...annotations, sampleAnnotation, ...xannotations],
  };
  // console.log("layout");
  // console.log(layout);

  //var config = { responsive: true };
  return { traces, layout };
}
