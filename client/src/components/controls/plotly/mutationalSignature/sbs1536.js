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

  const chunks = (a, size) =>
    Array.from(new Array(Math.ceil(a.length / size)), (_, i) =>
      a.slice(i * size, i * size + size)
    );
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

  const maxValMutation = Math.max(...data.map((o) => o.Mutations));
  console.log("maxValMutation:---");
  console.log(maxValMutation);

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

  const groupByMutationFront = data.reduce((groups, e, i) => {
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

  console.log("groupByMutationFront:---");
  console.log(groupByMutationFront);

  const mutationSumFront = Object.entries(groupByMutationFront).map(
    ([key, value]) => ({
      mutationType: key,
      contribution: value.reduce((a, e) => a + parseInt(e.contribution), 0),
    })
  );

  console.log("mutationSumFront:---");
  console.log(mutationSumFront);

  const arrayMutationSumFront = Object.values(mutationSumFront).flat();

  console.log("arrayMutationSumFront:---");
  console.log(arrayMutationSumFront);

  const groupByMutation2 = arrayMutationSumFront.reduce((groups, e, i) => {
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

  console.log("groupByMutation2:---");
  console.log(groupByMutation2);

  const heatmapY3 = [];
  const heatmapZ3 = [];
  const heatmapX3 = [];
  Object.entries(groupByMutation2).forEach(
    ([key, value], groupIndex, array) => {
      heatmapY3.push(Object.entries(value).map(([k, v]) => v.mutationType));
      heatmapZ3.push(
        Object.entries(value).map(([k, v]) => v.contribution / totalMutations)
      );
      heatmapX3.push(
        value.map(
          (e, i) =>
            array
              .slice(0, groupIndex)
              .reduce((x0, [_, sigs]) => x0 + sigs.length, 0) + i
        )
      );
    }
  );

  console.log("heatmapY3:---");
  console.log(heatmapY3);
  console.log("heatmapZ3:---");
  console.log(heatmapZ3);
  console.log("heatmapX3:---");
  console.log(heatmapX3);

  let heatmapY3_c = [
    heatmapY3[0][0].charAt(0) + "-- N",
    heatmapY3[0][16].charAt(0) + "-- N",
    heatmapY3[0][32].charAt(0) + "-- N",
    heatmapY3[0][48].charAt(0) + "-- N",
  ];
  console.log("heatmapY3_c:---");
  console.log(heatmapY3_c);
  let heatMapZ3_0 = chunks(heatmapZ3[0], 16);
  let heatMapZ3_1 = chunks(heatmapZ3[1], 16);
  let heatMapZ3_2 = chunks(heatmapZ3[2], 16);
  let heatMapZ3_3 = chunks(heatmapZ3[3], 16);
  let heatMapZ3_4 = chunks(heatmapZ3[4], 16);
  let heatMapZ3_5 = chunks(heatmapZ3[5], 16);

  heatmapZ3.forEach((item, index) => {});
  console.log(heatMapZ3_0);
  console.log(heatMapZ3_1);
  console.log(heatMapZ3_2);
  console.log(heatMapZ3_3);
  console.log(heatMapZ3_4);
  console.log(heatMapZ3_4);

  const heatMapZFinal3 = [
    heatMapZ3_0,
    heatMapZ3_1,
    heatMapZ3_2,
    heatMapZ3_3,
    heatMapZ3_4,
    heatMapZ3_5,
  ];

  const maxZ3 = Math.max(...heatMapZFinal3.flat(Infinity));
  console.log(maxZ3);

  const groupByMutationBack = data.reduce((groups, e, i) => {
    const mutation = e.MutationType.substring(1, e.MutationType.length);
    const signature = {
      mutationType: e.MutationType,
      contribution: e.Mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  console.log("groupByMutationBack:---");
  console.log(groupByMutationBack);

  const mutationSumBack = Object.entries(groupByMutationBack).map(
    ([key, value]) => ({
      mutationType: key,
      contribution: value.reduce((a, e) => a + parseInt(e.contribution), 0),
    })
  );

  console.log("mutationSumBack:---");
  console.log(mutationSumBack);

  const arrayMutationSumBack = Object.values(mutationSumBack).flat();

  console.log("arrayMutationSumBack:---");
  console.log(arrayMutationSumBack);

  const groupByMutation3 = arrayMutationSumBack.reduce((groups, e, i) => {
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

  console.log("groupByMutation3:---");
  console.log(groupByMutation3);

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

  const maxVal = Math.max(...flatSorted.map((o) => o.contribution));
  console.log("MaxVal");
  console.log(maxVal);

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

  console.log("heatmapZ");
  console.log(heatmapZ);

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
  // console.log(heatMapZ0);
  // console.log(heatMapZ1);
  // console.log(heatMapZ2);
  // console.log(heatMapZ3);
  // console.log(heatMapZ4);
  // console.log(heatMapZ5);

  // let heatMapZ0 = chunks(heatmapZ[0], 16);
  // let heatMapZ1 = chunks(heatmapZ[1], 16);
  // let heatMapZ2 = chunks(heatmapZ[2], 16);
  // let heatMapZ3 = chunks(heatmapZ[3], 16);
  // let heatMapZ4 = chunks(heatmapZ[4], 16);
  // let heatMapZ5 = chunks(heatmapZ[5], 16);

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

  const maxZ = Math.max(...heatMapZFinal.flat(Infinity));
  console.log(maxZ);

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

  const traceHeatMap = heatMapZFinal.map((num, index, array) => ({
    colorbar: { len: 0.35, y: 0.15, autotick: true, tick0: 0 },
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
    zmin: 0,
    zmax: maxZ,
    z: num,
    y: heatmapY,
    type: "heatmap",
    hoverongaps: false,
    xaxis: "x",
    yaxis: "y3",
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

  const traceHeatMap2 = heatMapZFinal3.map((num, index, array) => ({
    colorbar: { len: 0.35, y: 0.5, autotick: true, tick0: 0 },
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
    zmin: 0,
    zmax: maxZ3,
    z: num,
    y: heatmapY3_c,
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

  const traces = [...tracesBar, ...traceHeatMap, ...traceHeatMap2];

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
      rows: 3,
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
      anchor: "x",
      dtick: 1,
      tickfont: {
        size: 10,
      },
    },
    yaxis3: {
      autorange: true,
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
