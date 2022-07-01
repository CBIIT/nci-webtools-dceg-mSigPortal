export default function SBS192(data, sample) {
  const colors = {
    "C>A": "#03BCEE",
    "C>G": "black",
    "C>T": "#E32926",
    "T>A": "#CAC9C9",
    "T>C": "#A1CE63",
    "T>G": "#EBC6C4",
  };

  //   console.log("data--:");
  //   console.log(data);

  const arrayDataT = [];
  const arrayDataU = [];

  Object.values(data).forEach((group) => {
    if (group.mutationType.substring(0, 1) === "T") {
      arrayDataT.push(group);
    } else {
      arrayDataU.push(group);
    }
  });

  const numberWithCommas = (x) =>
    x.toString().replace(/\B(?<!\.\d*)(?=(\d{3})+(?!\d))/g, ",");
  const maxVal = Math.max(...data.map((o) => o.mutations));

  const totalMutations = data.reduce((a, e) => a + parseInt(e.mutations), 0);
  // group data by dominant mutation
  const groupByMutationT = arrayDataT.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.mutationType.match(mutationRegex)[1];
    const signature = {
      mutationType: e.mutationType.substring(2, e.mutationType.length),
      contribution: e.mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];

    return groups;
  }, {});
  const groupByMutationU = arrayDataU.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.mutationType.match(mutationRegex)[1];
    const signature = {
      mutationType: e.mutationType.substring(2, e.mutationType.length),
      contribution: e.mutations,
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
    marker: { color: "#004765" },
    x: flatSortedT.map((element, index, array) => index),
    y: flatSortedT.map((element, index, array) => element.contribution),

    hoverinfo: "x+y",
    showlegend: true,
  };

  const tracesU = {
    name: "Untranscribed Strand",
    type: "bar",
    marker: { color: "#E32925" },
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
  const xannotations = flatSortedU.map((num, index) => ({
    xref: "x",
    yref: "paper",
    xanchor: "bottom",
    yanchor: "bottom",
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
    align: "center",
    num: num,
    index: index,
    textangle: -90,
  }));

  const sampleAnnotation = {
    xref: "paper",
    yref: "paper",
    xanchor: "bottom",
    yanchor: "bottom",
    x: 0,
    y: 0.9,
    text:
      "<b>" +
      sample +
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
    hoverlabel: { bgcolor: "#FFF" },
    showlegend: true,
    legend: {
      x: 1,
      xanchor: "right",
      y: 1,
    },
    xaxis: {
      showticklabels: false,
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
      autorange: false,
      range: [0, maxVal + maxVal * 0.35],
      linecolor: "black",
      linewidth: 2,
      mirror: true,
    },

    shapes: [...shapes1, ...shapes2],
    annotations: [...annotations, sampleAnnotation, ...xannotations],
  };

  return { traces, layout };
}
