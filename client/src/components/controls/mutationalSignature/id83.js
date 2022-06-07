export default function ID83(data) {
  const colors = {
    "1:Del:C": "#FBBD6F",
    "1:Del:T": "#FE8002",
    "1:Ins:C": "#AEDD8A",
    "1:Ins:T": "#35A12E",
    "2:Del:R": "#FCC9B4",
    "3:Del:R": "#FB8969",
    "4:Del:R": "#F04432",
    "5:Del:R": "#BB1A1A",
    "2:Ins:R": "#CFDFF0",
    "3:Ins:R": "#93C3DE",
    "4:Ins:R": "#4B97C7",
    "5:Ins:R": "#1863AA",
    "2:Del:M": "#E1E1EE",
    "3:Del:M": "#B5B5D6",
    "4:Del:M": "#8482BC",
    "5:Del:M": "#62409A",
  };
  const annotationColors = {
    "1:Del:C": "black",
    "1:Del:T": "white",
    "1:Ins:C": "black",
    "1:Ins:T": "white",
    "2:Del:R": "black",
    "3:Del:R": "black",
    "4:Del:R": "black",
    "5:Del:R": "white",
    "2:Ins:R": "black",
    "3:Ins:R": "black",
    "4:Ins:R": "black",
    "5:Ins:R": "white",
    "2:Del:M": "blacl",
    "3:Del:M": "black",
    "4:Del:M": "black",
    "5:Del:M": "white",
  };
  let groupByFirstGroup = {},
    groupByMutationID = {},
    groupR = {},
    groupRDel = {},
    groupRIns = {},
    groupM = {},
    annotationIDTop = {},
    annotationIDBot = {},
    annotationsIDTopLabel = {},
    annotationsIDBotLabel = {},
    arrayID = [],
    arrayIDAnnotation = [],
    arrayIDAnnotationTop = [],
    arrayIDAnnXTop = [
      "1bp Deletion",
      "1bp Insertion",
      ">1bp Deletion at Repeats<br>(Deletion Length)",
      ">1bp Insertions at Repeats<br> (Insertion Length)",
      "Microhomology<br>(Deletion Length)",
    ],
    arrayIDAnnXBot = [
      "Homopolymer Length",
      "Homopolymer Length",
      "Number of Repeat Units",
      "Number of Repeat Units",
      "Microhimology Length",
    ],
    arrayIDAnnXLabel = [
      "1:Del:C:5",
      "1:Ins:C:5",
      "3:Del:R:5",
      "3:Ins:R:5",
      "4:Del:M:1",
    ];
  // group data by dominant mutation
  const groupByMutation = data.reduce((groups, e, i) => {
    const mutationRegex = /^.{0,7}/;
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

  console.log(groupByMutation);

  groupByFirstGroup = Object.fromEntries(
    Object.entries(groupByMutation).slice(0, 4)
  );

  groupByMutationID = data.reduce((groups, e) => {
    let mutationID;
    mutationID = e.MutationType.match(
      e.MutationType.substring(
        e.MutationType.length - 3,
        e.MutationType.length - 2
      )
    );
    const signature = {
      mutationType: e.MutationType,
      contribution: e.Contribution,
    };
    groups[mutationID] = groups[mutationID]
      ? [...groups[mutationID], signature]
      : [signature];
    return groups;
  }, {});

  groupR = groupByMutationID["R"].reduce((r, a) => {
    let m;
    m = a.mutationType.match(a.mutationType.substr(2, 3));
    const s = {
      mutationType: a.mutationType,
      contribution: a.contribution,
    };
    r[m] = r[m] ? [...r[m], a] : [s];
    return r;
  }, {});

  groupRDel = groupR["Del"].reduce((r, a) => {
    let m;
    m = a.mutationType.match(a.mutationType.substr(0, 7));
    const s = {
      mutationType: a.mutationType,
      contribution: a.contribution,
    };
    r[m] = r[m] ? [...r[m], a] : [s];
    return r;
  }, {});

  groupRIns = groupR["Ins"].reduce((r, a) => {
    let m;
    m = a.mutationType.match(a.mutationType.substr(0, 7));
    const s = {
      mutationType: a.mutationType,
      contribution: a.contribution,
    };
    r[m] = r[m] ? [...r[m], a] : [s];
    return r;
  }, {});

  groupM = groupByMutationID["M"].reduce((r, a) => {
    let m;
    m = a.mutationType.match(a.mutationType.substr(0, 7));
    const s = {
      mutationType: a.mutationType,
      contribution: a.contribution,
    };
    r[m] = r[m] ? [...r[m], a] : [s];
    return r;
  }, {});
  const arrayID1 = Object.keys(groupByFirstGroup).map(function (key) {
    return groupByFirstGroup[key];
  });
  const arrayID2 = Object.keys(groupRDel).map(function (key) {
    return groupRDel[key];
  });
  const arrayID3 = Object.keys(groupRIns).map(function (key) {
    return groupRIns[key];
  });
  const arrayID4 = Object.keys(groupM).map(function (key) {
    return groupM[key];
  });

  arrayID = [...arrayID1, ...arrayID2, ...arrayID3, ...arrayID4];
  console.log(arrayID);

  const flatSorted = Object.values(arrayID).flat();
  console.log(flatSorted);

  Object.values(arrayID).forEach((group) => {
    if (group.length > 1) {
      arrayIDAnnotationTop.push(
        group[Math.floor(group.length / 2)].mutationType
      );
    } else {
      arrayIDAnnotationTop.push(group[0].mutationType);
    }
    group.forEach((e) => {
      arrayIDAnnotation.push(e.mutationType);
    });
  });

  console.log(arrayIDAnnotation);
  console.log(arrayIDAnnotationTop);

  const traces = Object.entries(arrayID).map(
    ([mutation, signatures], groupIndex, array) => ({
      name: mutation,
      type: "bar",
      marker: {
        color:
          colors[
            signatures[0].mutationType.substring(
              0,
              signatures[0].mutationType.length - 2
            )
          ],
      },
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
    })
  );

  console.log(traces);

  const annotations1 = Object.entries(arrayID).map(
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
      y: 1.0,
      text:
        groupIndex < 4
          ? `<b>${signatures[0].mutationType.substring(
              signatures[0].mutationType.length - 3,
              signatures[0].mutationType.length - 2
            )}</b>`
          : `<b>${signatures[0].mutationType.substring(0, 1)}</b>`,
      showarrow: false,
      font: {
        size: 14,
        color:
          annotationColors[
            signatures[0].mutationType.substring(
              0,
              signatures[0].mutationType.length - 2
            )
          ],
      },
      align: "center",
      signatures: signatures,
      mutation: mutation,
      groupIndex: groupIndex,
    })
  );
  console.log(annotations1);

  const shapes1 = Object.entries(arrayID).map(
    ([mutation, signatures], groupIndex, array) => ({
      type: "rect",
      xref: "x",
      yref: "paper",
      x0: array
        .slice(0, groupIndex)
        .reduce((x0, [_, sigs]) => x0 + sigs.length, -0.4),
      x1: array
        .slice(0, groupIndex + 1)
        .reduce((x0, [_, sigs]) => x0 + sigs.length, -0.6),
      y0: 1.06,
      y1: 1,
      fillcolor:
        colors[
          signatures[0].mutationType.substring(
            0,
            signatures[0].mutationType.length - 2
          )
        ],
      line: {
        width: 0,
      },
      mutation: mutation,
      signature: signatures[0].mutationType.substring(
        0,
        signatures[0].mutationType.length - 2
      ),
    })
  );
  console.log(shapes1);

  const shapes2 = Object.entries(arrayID).map(
    ([mutation, signatures], groupIndex, array) => ({
      type: "rect",
      xref: "x",
      yref: "paper",
      x0: array
        .slice(0, groupIndex)
        .reduce((x0, [_, sigs]) => x0 + sigs.length, -0.4),
      x1: array
        .slice(0, groupIndex + 1)
        .reduce((x0, [_, sigs]) => x0 + sigs.length, -0.6),
      y0: -0.01,
      y1: -0.05,
      fillcolor:
        colors[
          signatures[0].mutationType.substring(
            0,
            signatures[0].mutationType.length - 2
          )
        ],
      line: {
        width: 0,
      },
      mutation: mutation,
      signature: signatures[0].mutationType.substring(
        0,
        signatures[0].mutationType.length - 2
      ),
    })
  );

  const layout = {
    xaxis: {
      title: "",
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
      title: "Number of Idels",
      autorange: true,
    },

    shapes: [...shapes1, ...shapes2],
    annotations: [...annotations1],
  };

  console.log(layout);
  return { traces, layout };
}
