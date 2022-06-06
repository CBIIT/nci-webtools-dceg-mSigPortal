export default function SBS96(data) {
  const colors = {
    "C>A": "#03BCEE",
    "C>G": "black",
    "C>T": "#E32926",
    "T>A": "#CAC9C9",
    "T>C": "#A1CE63",
    "T>G": "#EBC6C4",
  };

  // group data by dominant mutation
  const groupByMutation = data.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.MutationType.match(mutationRegex)[1];
    const signature = {
      mutationType: e.MutationType,
      contribution: e.Contribution,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});
  const flatSorted = Object.values(groupByMutation).flat();

  console.log(groupByMutation);
  console.log(flatSorted);

  const traces = Object.entries(groupByMutation).map(
    ([mutation, signatures], groupIndex) => ({
      name: mutation,
      type: "bar",
      marker: { color: colors[mutation] },
      //   x: signatures.map((e) => e.mutationType),
      x: signatures.map((e, i) => groupIndex * signatures.length + i),
      y: signatures.map((e) => e.contribution),
      hoverinfo: "x+y",
      showlegend: false,
    })
  );

  const annotations = Object.entries(groupByMutation).map(
    ([mutation, signatures], groupIndex) => ({
      xref: "x",
      yref: "paper",
      xanchor: "bottom",
      yanchor: "bottom",
      x: groupIndex * signatures.length + (signatures.length - 1) * 0.5,
      y: 1.04,
      text: `<b>${mutation}</b>`,
      showarrow: false,
      font: {
        size: 18,
      },
      align: "center",
    })
  );

  const shapes = Object.entries(groupByMutation).map(
    ([mutation, signatures], groupIndex) => ({
      type: "rect",
      xref: "x",
      yref: "paper",
      x0: groupIndex * signatures.length - 0.4,
      x1: groupIndex * signatures.length + signatures.length - 0.6,
      y0: 1.03,
      y1: 1,
      fillcolor: colors[mutation],
      line: {
        width: 0,
      },
    })
  );

  const layout = {
    xaxis: {
      title: "Substitution",
      showline: true,
      tickangle: -90,
      tickfont: {
        size: 10,
      },
      tickmode: "array",
      tickvals: flatSorted.map((_, i) => i),
      ticktext: flatSorted.map((e) => e.mutationType),
    },
    yaxis: {
      title: "Mutation Probability",
      autorange: true,
    },

    shapes: shapes,
    annotations: annotations,
  };

  return { traces, layout };
}
