export default function ID83(data, sample) {
  // console.log("data");
  // console.log(data);
  const colors = {
    '1:Del:C': '#FBBD6F',
    '1:Del:T': '#FE8002',
    '1:Ins:C': '#AEDD8A',
    '1:Ins:T': '#35A12E',
    '2:Del:R': '#FCC9B4',
    '3:Del:R': '#FB8969',
    '4:Del:R': '#F04432',
    '5:Del:R': '#BB1A1A',
    '2:Ins:R': '#CFDFF0',
    '3:Ins:R': '#93C3DE',
    '4:Ins:R': '#4B97C7',
    '5:Ins:R': '#1863AA',
    '2:Del:M': '#E1E1EE',
    '3:Del:M': '#B5B5D6',
    '4:Del:M': '#8482BC',
    '5:Del:M': '#62409A',
  };
  const annotationColors = {
    '1:Del:C': 'black',
    '1:Del:T': 'white',
    '1:Ins:C': 'black',
    '1:Ins:T': 'white',
    '2:Del:R': 'black',
    '3:Del:R': 'black',
    '4:Del:R': 'black',
    '5:Del:R': 'white',
    '2:Ins:R': 'black',
    '3:Ins:R': 'black',
    '4:Ins:R': 'black',
    '5:Ins:R': 'white',
    '2:Del:M': 'blacl',
    '3:Del:M': 'black',
    '4:Del:M': 'black',
    '5:Del:M': 'white',
  };

  const arrayIDAnnXTop = [
      '1bp Deletion',
      '1bp Insertion',
      '>1bp Deletion at Repeats<br>(Deletion Length)',
      '>1bp Insertions at Repeats<br> (Insertion Length)',
      'Microhomology<br>(Deletion Length)',
    ],
    arrayIDAnnXBot = [
      'Homopolymer Length',
      'Homopolymer Length',
      'Number of Repeat Units',
      'Number of Repeat Units',
      'Microhimology Length',
    ],
    arrayIDAnnXLabel = [5, 18.5, 35, 60, 76],
    arrayIDAnnotationTop = [],
    arrayIDAnnotationBot = [];

  const totalMutations = data.reduce((a, e) => a + parseInt(e.mutations), 0);
  const numberWithCommas = (x) =>
    x.toString().replace(/\B(?<!\.\d*)(?=(\d{3})+(?!\d))/g, ',');
  const maxVal = Math.max(...data.map((o) => o.mutations));

  // group data by dominant mutation
  const groupByMutation = data.reduce((groups, e, i) => {
    const mutationRegex = /^.{0,7}/;
    const mutation = e.mutationType.match(mutationRegex)[0];
    const signature = {
      mutationType: e.mutationType,
      contribution: e.mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  const groupByFirstGroup = Object.fromEntries(
    Object.entries(groupByMutation).slice(0, 4)
  );

  const groupByMutationID = data.reduce((groups, e) => {
    let mutationID;
    mutationID = e.mutationType.match(
      e.mutationType.substring(
        e.mutationType.length - 3,
        e.mutationType.length - 2
      )
    );
    const signature = {
      mutationType: e.mutationType,
      contribution: e.mutations,
    };
    groups[mutationID] = groups[mutationID]
      ? [...groups[mutationID], signature]
      : [signature];
    return groups;
  }, {});
  // console.log("groupByMutationID");
  // console.log(groupByMutationID);

  const groupR = groupByMutationID['R'].reduce((r, a) => {
    let m;
    m = a.mutationType.match(a.mutationType.substr(2, 3));
    const s = {
      mutationType: a.mutationType,
      contribution: a.contribution,
    };
    r[m] = r[m] ? [...r[m], a] : [s];
    return r;
  }, {});

  // console.log("groupR");
  // console.log(groupR);
  const groupRDel = groupR['Del'].reduce((r, a) => {
    let m;
    m = a.mutationType.match(a.mutationType.substr(0, 7));
    const s = {
      mutationType: a.mutationType,
      contribution: a.contribution,
    };
    r[m] = r[m] ? [...r[m], a] : [s];
    return r;
  }, {});

  const groupRIns = groupR['Ins'].reduce((r, a) => {
    let m;
    m = a.mutationType.match(a.mutationType.substr(0, 7));
    const s = {
      mutationType: a.mutationType,
      contribution: a.contribution,
    };
    r[m] = r[m] ? [...r[m], a] : [s];
    return r;
  }, {});

  const groupM = groupByMutationID['M'].reduce((r, a) => {
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

  const arrayID = [...arrayID1, ...arrayID2, ...arrayID3, ...arrayID4];

  const flatSorted = Object.values(arrayID).flat();

  Object.values(arrayID).forEach((group) => {
    if (group.length > 1) {
      arrayIDAnnotationTop.push(
        group[Math.floor(group.length / 2)].mutationType
      );
    } else {
      arrayIDAnnotationTop.push(group[0].mutationType);
    }
    group.forEach((e) => {
      arrayIDAnnotationBot.push(e.mutationType);
    });
  });

  const traces = Object.entries(arrayID).map(
    ([mutation, signatures], groupIndex, array) => ({
      name: mutation,
      type: 'bar',
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
      //text: signatures.map((e, i) => e.mutationType),
      //hovertemplate: "%{signatures.map((e, i) => e.mutationType)}, %{y}",
      hoverinfo: 'x+y',
      showlegend: false,
    })
  );

  const annotations1 = Object.entries(arrayID).map(
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
      y: 1.01,
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
      align: 'center',
      signatures: signatures,
      mutation: mutation,
      groupIndex: groupIndex,
    })
  );

  const annotations2 = arrayIDAnnotationBot.map((num, index) => ({
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: index,
    y: -0.1,
    text: '<b>' + num.substring(num.length - 1, num.length) + '</b>',
    showarrow: false,
    font: {
      size: 12,
    },
    align: 'center',
    num: num,
    index: index,
  }));

  const annotationsIDTopLabel = arrayIDAnnXLabel.map((num, index) => ({
    xref: 'x',
    yref: 'paper',
    x: num,
    xanchor: 'bottom',
    y: 1.07,
    yanchor: 'bottom',
    text: '<b>' + arrayIDAnnXTop[index] + '</b>',
    showarrow: false,
    font: {
      size: 14,
    },
    align: 'center',
  }));

  const annotationsIDBotLabel = arrayIDAnnXLabel.map((num, index) => ({
    xref: 'x',
    yref: 'paper',
    x: num,
    xanchor: 'bottom',
    y: -0.15,
    yanchor: 'bottom',
    text: '<b>' + arrayIDAnnXBot[index] + '</b>',
    showarrow: false,
    font: {
      size: 14,
    },
    align: 'center',
  }));

  const sampleAnnotation = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: 0,
    y: 0.92,
    text:
      '<b>' + sample + ': ' + numberWithCommas(totalMutations) + ' indels</b>',
    showarrow: false,
    font: {
      size: 18,
    },
    align: 'center',
  };

  const shapes1 = Object.entries(arrayID).map(
    ([mutation, signatures], groupIndex, array) => ({
      type: 'rect',
      xref: 'x',
      yref: 'paper',
      x0: array
        .slice(0, groupIndex)
        .reduce((x0, [_, sigs]) => x0 + sigs.length, -0.4),
      x1: array
        .slice(0, groupIndex + 1)
        .reduce((x0, [_, sigs]) => x0 + sigs.length, -0.6),
      y0: 1.07,
      y1: 1.01,
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

  const shapes2 = Object.entries(arrayID).map(
    ([mutation, signatures], groupIndex, array) => ({
      type: 'rect',
      xref: 'x',
      yref: 'paper',
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
    hoverlabel: { bgcolor: '#FFF' },
    height: 500,
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
      ticktext: flatSorted.map((e) => e.mutationType),
      linecolor: 'black',
      linewidth: 1,
      mirror: true,
    },
    yaxis: {
      title: 'Number of Idels',
      autorange: true,
      range: [0, maxVal + maxVal * 0.15],
      linecolor: 'black',
      linewidth: 1,
      mirror: true,
    },

    shapes: [...shapes1, ...shapes2],
    annotations: [
      ...annotations1,
      ...annotations2,
      ...annotationsIDTopLabel,
      ...annotationsIDBotLabel,
      sampleAnnotation,
    ],
  };

  var config = { responsive: true };
  //console.log("layout");
  //console.log(layout);
  return { traces, layout, config };
}
