export default function SBS384(data, sample) {
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
  const numberWithCommas = (x) =>
    x.toString().replace(/\B(?<!\.\d*)(?=(\d{3})+(?!\d))/g, ",");

  const groupByMutation = data.reduce((groups, e, i) => {
    const mutation = e.mutationType.substring(2, e.mutationType.length);
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

  console.log("groupByMutation:");
  console.log(groupByMutation);
  console.log("FlatSorted");
  console.log(flatSorted);
  //group data by 1st letter

  const dataT = [];
  const dataU = [];

  Object.entries(flatSorted).forEach(([key, value], groupIndex, array) => {
    if (value.mutationType.substring(0, 1) === "T") {
      dataT.push(value);
    } else if (value.mutationType.substring(0, 1) === "U") {
      dataU.push(value);
    }
  });

  const totalMutations = [...dataT, ...dataU].reduce(
    (a, e) => a + parseInt(e.contribution),
    0
  );

  const dataUT = [...dataT, ...dataU];
  const maxVal = Math.max(...dataUT.map((o) => o.contribution));
  console.log("maxVal--:");
  console.log(maxVal);

  console.log("dataT");
  console.log(dataT);
  console.log("dataU");
  console.log(dataU);
  const groupByMutationU = dataU.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.mutationType.match(mutationRegex)[1];
    const signature = {
      mutationType: e.mutationType.substring(2, e.mutationType.length),
      contribution: e.contribution,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});
  console.log("groupByMutationU");
  console.log(groupByMutationU);
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
  console.log("flatSortedU");
  console.log(flatSortedU);

  const groupByMutationT = dataT.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.mutationType.match(mutationRegex)[1];
    const signature = {
      mutationType: e.mutationType.substring(2, e.mutationType.length),
      contribution: e.contribution,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  console.log("groupByMutationT");
  console.log(groupByMutationT);
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
  console.log("flatSortedT");
  console.log(flatSortedT);

  console.log("dataT");
  console.log(dataT);
  console.log("dataU");
  console.log(dataU);

  const tracesT = {
    name: "Transcrribed Strand",
    type: "bar",
    marker: { color: "#004765" },
    x: flatSortedT.map((element, index, array) => index),
    y: flatSortedT.map((element, index, array) => element.contribution),

    hoverinfo: "x+y",
    showlegend: true,
  };

  console.log(tracesT);
  const tracesU = {
    name: "Untranscribed",
    type: "bar",
    marker: { color: "#E32925" },
    x: flatSortedU.map((element, index, array) => index),
    y: flatSortedU.map((element, index, array) => element.contribution),

    hoverinfo: "x+y",
    showlegend: true,
  };
  console.log(tracesU);

  const traces = [tracesT, tracesU];

  const annotations = Object.entries(groupByMutationT).map(
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

  console.log("annotations:");
  console.log(annotations);

  const xannotations = flatSortedT.map((num, index) => ({
    xref: "x",
    yref: "paper",
    xanchor: "bottom",
    yanchor: "bottom",
    x: index,
    y: -0.065,
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

  console.log("xannotations");
  console.log(xannotations);

  const sampleAnnotation = {
    xref: "paper",
    yref: "paper",
    xanchor: "bottom",
    yanchor: "bottom",
    x: 0,
    y: 0.92,
    text:
      "<b>" + sample + ": " + numberWithCommas(totalMutations) + " subs </b>",
    showarrow: false,
    font: {
      size: 18,
    },
    align: "center",
  };

  const shapes1 = Object.entries(groupByMutationT).map(
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
      mutation: mutation,
    })
  );
  const shapes2 = Object.entries(groupByMutationT).map(
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

  const layout = {
    hoverlabel: { bgcolor: "#FFF" },
    bargap: 0.3,
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
      tickvals: [...flatSortedT.map((_, i) => i)],
      ticktext: [...flatSortedT.map((e) => e.mutationType)],
      linecolor: "black",
      linewidth: 2,
      mirror: true,
    },
    yaxis: {
      title: "Number of Single Base Substitutions",
      autorange: false,
      range: [0, maxVal + maxVal * 0.15],
      linecolor: "black",
      linewidth: 2,
      mirror: true,
      categoryorder: "category descending",
    },
    shapes: [...shapes1, ...shapes2],
    annotations: [...annotations, sampleAnnotation, ...xannotations],
  };
  // console.log("layout");
  // console.log(layout);

  return { traces, layout };
}
