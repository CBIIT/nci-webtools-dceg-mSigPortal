export default function DBS78(data) {
  const colors = {
    "AC>": "#09BCED",
    "AT>": "#0266CA",
    "CC>": "#9FCE62",
    "CG>": "#006501",
    "CT>": "#FF9898",
    "GC>": "#E22925",
    "TA>": "#FEB065",
    "TC>": "#FD8000",
    "TG>": "#CB98FD",
    "TT>": "#4C0299",
  };

  // group data by dominant mutation
  const groupByMutation = data.reduce((groups, e, i) => {
    const mutationRegex = /^.{0,3}/;
    const mutation = e.MutationType.match(mutationRegex)[0];
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

  const xValues = flatSorted.map((e, i) => i);
  console.log(xValues);
  console.log(groupByMutation);
  console.log(flatSorted);
  const traces = Object.entries(groupByMutation).map(
    ([mutation, signatures], groupIndex, array) => ({
      name: mutation,
      type: "bar",
      marker: { color: colors[mutation] },
      //x: signatures.map((e) => e.mutationType),
      x: signatures.map(
        (e, i) =>
          array
            .slice(0, groupIndex)
            .reduce((x0, [_, sigs]) => x0 + sigs.length, 0) + i
      ),
      y: signatures.map((e) => e.contribution),
      hoverinfo: "x+y",
      showlegend: false,
    })
  );
  console.log(traces);
  const annotations = Object.entries(groupByMutation).map(
    ([mutation, signatures], groupIndex, array) => ({
      xref: "x",
      yref: "paper",
      xanchor: "bottom",
      yanchor: "bottom",
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
      align: "center",
    })
  );
  console.log(annotations);
  const shapes = Object.entries(groupByMutation).map(
    ([mutation, signatures], groupIndex, array) => ({
      type: "rect",
      xref: "x",
      yref: "paper",
      x0: array
        .slice(0, groupIndex)
        .reduce((x0, [_, signatures]) => x0 + signatures.length, -0.4),
      // x0: groupIndex * 16 - 0.4,
      x1: array
        .slice(0, groupIndex + 1)
        .reduce((x0, [_, signatures]) => x0 + signatures.length, -0.6),
      // x1: groupIndex * 16 + signatures.length - 0.6,
      y0: 1.03,
      y1: 1,
      fillcolor: colors[mutation],
      line: {
        width: 0,
      },
    })
  );
  console.log(shapes);

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
