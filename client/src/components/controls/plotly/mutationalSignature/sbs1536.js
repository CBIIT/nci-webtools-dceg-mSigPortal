export default function SBS96(data, sample) {
  const colors = {
    "C>A": "#03BCEE",
    "C>G": "black",
    "C>T": "#E32926",
    "T>A": "#CAC9C9",
    "T>C": "#A1CE63",
    "T>G": "#EBC6C4",
  };
  const numberWithCommas = (x) =>
    x.toString().replace(/\B(?<!\.\d*)(?=(\d{3})+(?!\d))/g, ",");

  const totalMutations = data.reduce((a, e) => a + parseInt(e.Mutations), 0);
  const maxVal = Math.max(...data.map((o) => o.Mutations));

  console.log("data--:");
  console.log(data);

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

  console.log("groupByMutation:---");
  console.log(groupByMutation);

  const groupByMutationInner = data.reduce((groups, e, i) => {
    const mutation = e.MutationType.substring(1, 8);

    const signature = {
      mutationType: e.MutationType,
      contribution: e.Mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  //   console.log("groupByMutationInner:---");
  //   console.log(groupByMutationInner);

  const groupByMutationOuter = data.reduce((groups, e, i) => {
    const mutation =
      e.MutationType[0] + e.MutationType[e.MutationType.length - 1];
    const signature = {
      mutationType: e.MutationType,
      contribution: e.Mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  console.log("groupByMutationOuter:---");
  console.log(groupByMutationOuter);

  const groupByMutationOuterInner = data.reduce((groups, e, i) => {
    const mutation =
      e.MutationType[0] +
      e.MutationType.substring(2, 7) +
      e.MutationType[e.MutationType.length - 1];
    const signature = {
      mutationType: e.MutationType,
      contribution: e.Mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  console.log("groupByMutationOuterInner:---");
  console.log(groupByMutationOuterInner);

  const groupByMutationTN = data.reduce((groups, e, i) => {
    const mutation = e.MutationType[0] + e.MutationType.substring(2, 7);
    const signature = {
      mutationType: e.MutationType,
      contribution: e.Mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  console.log("groupByMutationTN:---");
  console.log(groupByMutationTN);

  //   const testObj = Object.fromEntries(
  //     Object.entries(groupByMutationTN).map(([key, value]) => {
  //       const sumVal = value.reduce((mutationSum, arrayVal) => {
  //         const muType = arrayVal.mutationType;
  //         const currMutationName = muType.substring(0, muType.length - 1);
  //         console.log(mutationSum);
  //         console.log(arrayVal);
  //         console.log(arrayVal.contribution);
  //         //   if (mutationSum[currMutationName]) {
  //         //     mutationSum[currMutationName] =
  //         //       mutationSum[currMutationName] + arrayVal.contribution;
  //         //   } else {
  //         //     mutationSum[currMutationName] = arrayVal.contribution;
  //         //   }
  //         // }, {});
  //       });
  //       return [key, sumVal];
  //     })
  //   );

  //   console.log("testObj:---");
  //   console.log(testObj);

  const groupByMutationSum4 = data.reduce((groups, e, i) => {
    const mutation = e.MutationType.substring(0, e.MutationType.length - 1);
    const signature = {
      mutationType: e.MutationType,
      contribution: e.Mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  console.log("groupByMutationSum4:---");
  console.log(groupByMutationSum4);

  const groupByMutationA = data.reduce((groups, e, i) => {
    const mutation = e.MutationType.substring(0, 8);
    const signature = {
      mutationType: e.MutationType,
      contribution: e.Mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  console.log("groupByMutationA:---");
  console.log(groupByMutationA);

  const totalMutationsGroup = Object.entries(groupByMutationInner).map(
    ([mutation, signatures], groupIndex, array) => ({
      mutationType: mutation,
      signatures: signatures,
      groupIndex: groupIndex,
      array: array,
      total: signatures.reduce((a, e) => a + parseInt(e.contribution), 0),
    })
  );

  console.log("totalMutationsGroup");
  console.log(totalMutationsGroup);

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
  console.log("groupByTotal");
  console.log(groupByTotal);

  const flatSorted = Object.values(groupByTotal).flat();
  const mutationTitle = Object.keys(groupByTotal).flat();

  console.log("FlatSorted--");
  console.log(flatSorted);
  console.log("mutationTitle--");
  console.log(mutationTitle);

  const arrayDataX1 =
    groupByMutationOuter[Object.keys(groupByMutationOuter)[0]];
  console.log(arrayDataX1);

  const groupByMutation1 = arrayDataX1.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.mutationType.match(mutationRegex)[1];

    const signature = {
      mutationType: e.mutationType,
      contribution: e.contribution,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});
  //   console.log("groupByMutation1");
  //   console.log(groupByMutation1);

  const heatmapY = [];
  const heatmapZ = [];
  const heatmapX = [];
  Object.entries(groupByMutationOuter).forEach(
    ([key, value], groupIndex, array) => {
      heatmapY.push(key.charAt(0) + " -- " + key.charAt(key.length - 1));
      heatmapZ.push(
        Object.entries(value).map(([k, v]) => v.contribution / totalMutations)
      );
      heatmapX.push(
        value.map(
          (e, i) =>
            array
              .slice(0, groupIndex)
              .reduce((x0, [_, sigs]) => x0 + sigs.length, 0) + i
        )
      );
    }
  );

  console.log("heatmapY");
  //   console.log(heatmapY);
  console.log("heatmapZ");
  console.log(heatmapZ);
  //   console.log("heatmapX");
  //   console.log(heatmapX);

  //const copyHeatMapZ = [...heatmapZ];
  let heatMapZ0 = [];
  let heatMapZ1 = [];
  let heatMapZ2 = [];
  let heatMapZ3 = [];
  let heatMapZ4 = [];
  let heatMapZ5 = [];

  heatmapZ.forEach((item, index) => {
    // console.log(item);
    // console.log(index);

    heatMapZ0.push(item.slice().splice(0, 16));
    heatMapZ1.push(item.slice().splice(16, 16));
    heatMapZ2.push(item.slice().splice(32, 16));
    heatMapZ3.push(item.slice().splice(48, 16));
    heatMapZ4.push(item.slice().splice(64, 16));
    heatMapZ5.push(item.slice().splice(80, 16));
    //console.log(splitToChunks(item, 6));
  });
  console.log(heatMapZ0);
  console.log(heatMapZ1);
  console.log(heatMapZ2);
  console.log(heatMapZ3);
  console.log(heatMapZ4);
  console.log(heatMapZ5);

  const heatMapZFinal = [
    heatMapZ0,
    heatMapZ1,
    heatMapZ2,
    heatMapZ3,
    heatMapZ4,
    heatMapZ5,
    //mutationTitle,
  ];
  console.log("heatMapZFinal");
  console.log(heatMapZFinal);

  const tracesBar = Object.entries(groupByTotal).map(
    ([mutation, signatures], groupIndex, array) => ({
      name: mutation,
      type: "bar",
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
      hoverinfo: "x+y",
      showlegend: false,
      array: array,
      mutation: mutation,
      signatures: signatures,
    })
  );
  console.log("tracesBar:");
  console.log(tracesBar);

  const tracesHeat = Object.entries(groupByMutation).map(
    ([mutation, signatures], groupIndex, array) => ({
      name: mutation,
      type: "heatmap",
      marker: { color: colors[mutation] },
      //   x: signatures.map((e) => e.mutationType),
      //x: signatures.map((e, i) => groupIndex * signatures.length + i),
      x: signatures.map(
        (e, i) =>
          array
            .slice(0, groupIndex)
            .reduce((x0, [_, sigs]) => x0 + sigs.length, 0) + i
      ),
      y: heatmapY,
      z: Object.entries(signatures).map(([k, v]) => v.contribution),
      hoverinfo: "x+y",
      showlegend: false,
      array: array,
      mutation: mutation,
      signatures: signatures,
      xaxis: "x",
      yaxis: "y2",
      colorbar: { len: 0.5, y: 0.2 },
    })
  );

  //   const tracesHeat1 = [
  //     {
  //       colorbar: { len: 0.5, y: 0.2 },
  //       z: heatmapZ,
  //       x: heatmapX[0],
  //       y: heatmapY,
  //       type: "heatmap",
  //       hoverongaps: false,
  //       xaxis: "x",
  //       yaxis: "y2",
  //     },
  //   ];
  console.log("tracesHeat:");
  console.log(tracesHeat);

  //   console.log("tracesHeat1:");
  //   console.log(tracesHeat1);

  const traceHeatMap = heatMapZFinal.map((num, index, array) => ({
    colorbar: { len: 0.5, y: 0.2, autotick: true, tick0: 0, dtick: 0.005 },
    colorscale: [
      [0, "rgb(56,56,156"],
      [0.2, "rgb(56,56,156"],
      [0.2, "rgb(106,106,128"],
      [0.4, "rgb(106,106,128"],
      [0.4, "rgb(155,146,98"],
      [0.6, "rgb(155,146,98"],
      [0.6, "rgb(205,186,69"],
      [0.8, "rgb(205,186,69"],
      [0.8, "rgb(255,255,39)"],
      [1, "rgb(255,255,39)"],
    ],
    z: num,
    y: heatmapY,
    type: "heatmap",
    hoverongaps: false,
    xaxis: "x",
    yaxis: "y2",
    num: num,

    x: num.map(
      (e, i) =>
        array.slice(0, index).reduce((x0, [_, sigs]) => x0 + sigs.length, 0) + i
    ),
    hovertemplate:
      "x: %{x}<br>" + "y: %{y}<br>" + "Value: %{z}" + "<extra></extra>",
  }));

  console.log("traceHeatMap");
  console.log(traceHeatMap);

  const traces = [...tracesBar, ...traceHeatMap];

  console.log("traces:");
  console.log(traces);
  const annotations = Object.entries(groupByTotal).map(
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

  console.log("annotation:");
  console.log(annotations);

  const sampleAnnotation = {
    xref: "paper",
    yref: "paper",
    xanchor: "bottom",
    yanchor: "bottom",
    x: 0,
    y: 0.9,
    text:
      "<b>" + sample + ": " + numberWithCommas(totalMutations) + " subs </b>",
    showarrow: false,
    font: {
      size: 18,
    },
    align: "center",
  };

  const shapes = Object.entries(groupByTotal).map(
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

  const xannotations = flatSorted.map((num, index) => ({
    xref: "x",
    yref: "paper",
    xanchor: "bottom",
    yanchor: "bottom",
    x: index,
    y: -0.1,
    text: num.mutationType.replace(/\[(.*)\]/, "-"),
    showarrow: false,
    font: {
      size: 10,
      //   color: colors[num.mutationType.substring(2, 5)],
    },
    align: "center",
    num: num,
    index: index,
    textangle: -90,
  }));
  const layout = {
    hoverlabel: { bgcolor: "#FFF" },
    grid: {
      rows: 2,
      columns: 1,
    },
    xaxis: {
      showticklabels: false,
      showline: true,
      tickangle: -90,
      tickfont: {
        size: 10,
      },
      tickmode: "array",
      tickvals: flatSorted.map((_, i) => i),
      ticktext: flatSorted.map((e) => e.mutationType),
      linecolor: "black",
      linewidth: 1,
      mirror: true,
    },
    yaxis: {
      //   title: "Number of Single Base Substitutions",
      autorange: false,
      range: [0, maxVal + maxVal * 0.2],
      linecolor: "black",
      linewidth: 1,
      mirror: true,
    },
    yaxis2: {
      autorange: true,
      linecolor: "black",
      linewidth: 1,
      mirror: true,
      anchor: "x",
      dtick: 1,
      tickfont: {
        size: 10,
      },
    },

    shapes: shapes,
    annotations: [...annotations, sampleAnnotation, ...xannotations],
  };
  console.log("layout");
  console.log(layout);

  return { traces, layout };
}
