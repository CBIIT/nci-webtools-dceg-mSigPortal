export default function SBS192(data) {
  const colors = {
    "C>A": "#03BCEE",
    "C>G": "black",
    "C>T": "#E32926",
    "T>A": "#CAC9C9",
    "T>C": "#A1CE63",
    "T>G": "#EBC6C4",
  };
  function sum(arr) {
    const sum = arr.reduce((a, b) => a + b, 0);
    return sum || 0;
  }
  console.log("data--:");
  console.log(data);

  const arrayDataT = [];
  const arrayDataU = [];

  function dynamicSort(property) {
    var sortOrder = 1;

    if (property[0] === "-") {
      sortOrder = -1;
      property = property.substr(1);
    }

    return function (a, b) {
      if (sortOrder == -1) {
        return b[property].localeCompare(a[property]);
      } else {
        return a[property].localeCompare(b[property]);
      }
    };
  }

  Object.values(data).forEach((group) => {
    if (group.MutationType.substring(0, 1) === "T") {
      arrayDataT.push(group);
    } else {
      arrayDataU.push(group);
    }
  });

  console.log("group T");
  console.log(arrayDataT);
  console.log("group U");
  console.log(arrayDataU);
  // group data by dominant mutation
  const groupByMutation = data.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.MutationType.match(mutationRegex)[1];
    const signature = {
      mutationType: e.MutationType,
      contribution: e.Mutations,
    };

    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});
  const groupByMutationT = arrayDataT.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.MutationType.match(mutationRegex)[1];
    const signature = {
      mutationType: e.MutationType.substring(2, e.MutationType.length),
      contribution: e.Mutations,
    };
    console.log(signature);

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
    console.log(signature);
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

  const groupByMutationTS = flatSortedT.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.mutationType.match(mutationRegex)[1];
    const signature = {
      mutationType: e.mutationType,
      contribution: e.contribution,
    };
    console.log(signature);

    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];

    return groups;
  }, {});

  const groupByMutationUS = flatSortedU.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.mutationType.match(mutationRegex)[1];
    const signature = {
      mutationType: e.mutationType,
      contribution: e.contribution,
    };
    console.log(signature);

    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];

    return groups;
  }, {});
  console.log("groupByMutation T:");
  console.log(groupByMutationT);
  console.log("groupByMutation U:");
  console.log(groupByMutationU);
  console.log("FlatSorted T");
  console.log(flatSortedT);
  console.log("FlatSorted U");
  console.log(flatSortedU);
  const tracesT = Object.entries(groupByMutationTS).map(
    ([mutation, signatures], groupIndex, array) => ({
      name: mutation,
      type: "bar",
      marker: { color: "blue" },
      //   x: signatures.map((e) => e.mutationType),
      //x: signatures.map((e, i) => groupIndex * signatures.length + i),
      x: signatures.map(
        (e, i) =>
          array
            .slice(0, groupIndex)
            .reduce((x0, [_, sigs]) => x0 + sigs.length, 0) + i
      ),
      y: signatures.map((e) => e.contribution),
      sum: sum(signatures.map((e) => e.contribution)),
      hoverinfo: "x+y",
      showlegend: false,
      array: array,
    })
  );
  console.log("traces T:");
  console.log(tracesT);

  const tracesU = Object.entries(groupByMutationUS).map(
    ([mutation, signatures], groupIndex, array) => ({
      name: mutation,
      type: "bar",
      marker: { color: "red" },
      //   x: signatures.map((e) => e.mutationType),
      //x: signatures.map((e, i) => groupIndex * signatures.length + i),
      x: signatures.map(
        (e, i) =>
          array
            .slice(0, groupIndex)
            .reduce((x0, [_, sigs]) => x0 + sigs.length, 0) + i
      ),
      y: signatures.map((e) => e.contribution),
      sum: sum(signatures.map((e) => e.contribution)),
      hoverinfo: "x+y",
      showlegend: false,
      array: array,
    })
  );
  console.log("traces U:");
  console.log(tracesU);

  const annotations = Object.entries(groupByMutationUS).map(
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
    y: 0.85,
    text: "<b>" + data[0].Sample + "</b>",
    showarrow: false,
    font: {
      size: 18,
    },
    align: "center",
  };

  const shapes1 = Object.entries(groupByMutationUS).map(
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

  const shapes2 = Object.entries(groupByMutationUS).map(
    ([mutation, signatures], groupIndex, array) => ({
      type: "rect",
      xref: "x",
      yref: "paper",
      //x0: signatures.map((e) => e.mutationType)[0],
      y0: 1,
      //x1: signatures.map((e) => e.mutationType)[signatures.length - 1],
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
  const traces = [...tracesT, ...tracesU];
  console.log("traces:");
  console.log(traces);

  const layout = {
    barmode: "group",
    xaxis: {
      showticklabels: true,
      showline: true,
      tickangle: -90,
      tickfont: {
        size: 10,
      },
      tickmode: "array",
      tickvals: [
        ...flatSortedT.map((_, i) => i),
        ...flatSortedU.map((_, i) => i),
      ],
      ticktext: [
        ...flatSortedT.map((e) => e.mutationType),
        ...flatSortedU.map((e) => e.mutationType),
      ],
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
  console.log("layout");
  console.log(layout);

  return { traces, layout };
}
