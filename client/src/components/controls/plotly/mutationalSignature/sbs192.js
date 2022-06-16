export default function SBS192(data) {
  const colors = {
    "C>A": "#03BCEE",
    "C>G": "black",
    "C>T": "#E32926",
    "T>A": "#CAC9C9",
    "T>C": "#A1CE63",
    "T>G": "#EBC6C4",
  };

  console.log("data--:");
  console.log(data);

  const arrayDataT = [];
  const arrayDataU = [];

  Object.values(data).forEach((group) => {
    if (group.MutationType.substring(0, 1) === "T") {
      arrayDataT.push(group);
    } else {
      arrayDataU.push(group);
    }
  });

  const numberWithCommas = (x) =>
    x.toString().replace(/\B(?<!\.\d*)(?=(\d{3})+(?!\d))/g, ",");

  const totalMutations = data.reduce((a, e) => a + parseInt(e.Mutations), 0);
  // group data by dominant mutation
  const groupByMutationT = arrayDataT.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.MutationType.match(mutationRegex)[1];
    const signature = {
      mutationType: e.MutationType.substring(2, e.MutationType.length),
      contribution: e.Mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];

    return groups;
  }, {});
  const groupByMutationU = arrayDataU.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.MutationType.match(mutationRegex)[1];
    const signature = {
      mutationType: e.MutationType.substring(2, e.MutationType.length),
      contribution: e.Mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

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
    name: "Transcrribed Strand",
    type: "bar",
    marker: { color: "blue" },
    x: flatSortedT.map((element, index, array) => index),
    y: flatSortedT.map((element, index, array) => element.contribution),

    hoverinfo: "x+y",
    showlegend: true,
  };

  const tracesU = {
    name: "Untranscribed Strand",
    type: "bar",
    marker: { color: "red" },
    x: flatSortedU.map((element, index, array) => index),
    y: flatSortedU.map((element, index, array) => element.contribution),

    hoverinfo: "x+y",
    showlegend: true,
  };

  const annotations = Object.entries(groupByMutationU).map(
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
      text: `<b>${mutation}</b>`,
      showarrow: false,
      font: {
        size: 18,
      },
      align: "center",
    })
  );

  const sampleAnnotation = {
    xref: "paper",
    yref: "paper",
    xanchor: "bottom",
    yanchor: "bottom",
    x: 0,
    y: 0.9,
    text:
      "<b>" +
      data[0].Sample +
      ": " +
      numberWithCommas(totalMutations) +
      " transcribed subs</b>",
    showarrow: false,
    font: {
      size: 18,
    },
    align: "center",
  };

  const shapes1 = Object.entries(groupByMutationU).map(
    ([mutation, _], groupIndex, array) => ({
      type: "rect",
      xref: "x",
      yref: "paper",
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
      type: "rect",
      xref: "x",
      yref: "paper",
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
      opacity: 0.2,
    })
  );
  const traces = [tracesT, tracesU];

  const layout = {
    showlegend: true,
    legend: {
      x: 1,
      xanchor: "right",
      y: 1,
    },
    xaxis: {
      showticklabels: true,
      showline: true,
      tickangle: -90,
      tickfont: {
        size: 10,
      },
      tickmode: "array",
      tickvals: [...flatSortedU.map((_, i) => i)],
      ticktext: [...flatSortedU.map((e) => e.mutationType)],
      linecolor: "black",
      linewidth: 2,
      mirror: true,
    },
    yaxis: {
      title: "Number of Single Base Substitutions",
      autorange: true,
      linecolor: "black",
      linewidth: 2,
      mirror: true,
    },

    shapes: [...shapes1, ...shapes2],
    annotations: [...annotations, sampleAnnotation],
  };

  return { traces, layout };
}
