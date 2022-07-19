export default function ID83(data, sample) {
  const colors = {
    '1:Del:C': { shape: '#FBBD6F', text: 'black' },
    '1:Del:T': { shape: '#FE8002', text: 'white' },
    '1:Ins:C': { shape: '#AEDD8A', text: 'black' },
    '1:Ins:T': { shape: '#35A12E', text: 'white' },
    '2:Del:R': { shape: '#FCC9B4', text: 'black' },
    '3:Del:R': { shape: '#FB8969', text: 'black' },
    '4:Del:R': { shape: '#F04432', text: 'black' },
    '5:Del:R': { shape: '#BB1A1A', text: 'white' },
    '2:Ins:R': { shape: '#CFDFF0', text: 'black' },
    '3:Ins:R': { shape: '#93C3DE', text: 'black' },
    '4:Ins:R': { shape: '#4B97C7', text: 'black' },
    '5:Ins:R': { shape: '#1863AA', text: 'white' },
    '2:Del:M': { shape: '#E1E1EE', text: 'blacl' },
    '3:Del:M': { shape: '#B5B5D6', text: 'black' },
    '4:Del:M': { shape: '#8482BC', text: 'black' },
    '5:Del:M': { shape: '#62409A', text: 'white' },
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

  const totalMutations = data.reduce(
    (total, indel) =>
      total + indel.data.reduce((indelSum, e) => indelSum + e.mutations, 0),
    0
  );
  const maxMutation = Math.max(
    ...data.map((indel) => indel.data.map((e) => e.mutations)).flat()
  );
  const indelNames = data
    .map((indel) =>
      indel.data.map((e) => ({
        indel: indel.indel,
        index: e.mutationType.slice(-1),
      }))
    )
    .flat();

  // const groupByFirstGroup = data.slice(0, 4);

  // const groupByMutationID = data.reduce((groups, e) => {
  //   let mutationID;
  //   mutationID = e.mutationType.match(
  //     e.mutationType.substring(
  //       e.mutationType.length - 3,
  //       e.mutationType.length - 2
  //     )
  //   );
  //   const signature = {
  //     mutationType: e.mutationType,
  //     contribution: e.mutations,
  //   };
  //   groups[mutationID] = groups[mutationID]
  //     ? [...groups[mutationID], signature]
  //     : [signature];
  //   return groups;
  // }, {});
  // console.log("groupByMutationID");
  // console.log(groupByMutationID);

  // const groupR = groupByMutationID['R'].reduce((r, a) => {
  //   let m;
  //   m = a.mutationType.match(a.mutationType.substr(2, 3));
  //   const s = {
  //     mutationType: a.mutationType,
  //     contribution: a.contribution,
  //   };
  //   r[m] = r[m] ? [...r[m], a] : [s];
  //   return r;
  // }, {});

  // // console.log("groupR");
  // // console.log(groupR);
  // const groupRDel = groupR['Del'].reduce((r, a) => {
  //   let m;
  //   m = a.mutationType.match(a.mutationType.substr(0, 7));
  //   const s = {
  //     mutationType: a.mutationType,
  //     contribution: a.contribution,
  //   };
  //   r[m] = r[m] ? [...r[m], a] : [s];
  //   return r;
  // }, {});

  // const groupRIns = groupR['Ins'].reduce((r, a) => {
  //   let m;
  //   m = a.mutationType.match(a.mutationType.substr(0, 7));
  //   const s = {
  //     mutationType: a.mutationType,
  //     contribution: a.contribution,
  //   };
  //   r[m] = r[m] ? [...r[m], a] : [s];
  //   return r;
  // }, {});

  // const groupM = groupByMutationID['M'].reduce((r, a) => {
  //   let m;
  //   m = a.mutationType.match(a.mutationType.substr(0, 7));
  //   const s = {
  //     mutationType: a.mutationType,
  //     contribution: a.contribution,
  //   };
  //   r[m] = r[m] ? [...r[m], a] : [s];
  //   return r;
  // }, {});
  // const arrayID1 = Object.keys(groupByFirstGroup).map(function (key) {
  //   return groupByFirstGroup[key];
  // });
  // const arrayID2 = Object.keys(groupRDel).map(function (key) {
  //   return groupRDel[key];
  // });
  // const arrayID3 = Object.keys(groupRIns).map(function (key) {
  //   return groupRIns[key];
  // });
  // const arrayID4 = Object.keys(groupM).map(function (key) {
  //   return groupM[key];
  // });

  // const arrayID = [...arrayID1, ...arrayID2, ...arrayID3, ...arrayID4];

  // const flatSorted = Object.values(arrayID).flat();

  // Object.values(arrayID).forEach((group) => {
  //   if (group.length > 1) {
  //     arrayIDAnnotationTop.push(
  //       group[Math.floor(group.length / 2)].mutationType
  //     );
  //   } else {
  //     arrayIDAnnotationTop.push(group[0].mutationType);
  //   }
  //   group.forEach((e) => {
  //     arrayIDAnnotationBot.push(e.mutationType);
  //   });
  // });

  // const traces = Object.entries(arrayID).map(
  //   ([mutation, signatures], groupIndex, array) => ({
  //     name: mutation,
  //     type: 'bar',
  //     marker: {
  //       color:
  //         colors[
  //           signatures[0].mutationType.substring(
  //             0,
  //             signatures[0].mutationType.length - 2
  //           )
  //         ],
  //     },
  //     //   x: signatures.map((e) => e.mutationType),
  //     //x: signatures.map((e, i) => groupIndex * signatures.length + i),
  //     x: signatures.map(
  //       (e, i) =>
  //         array
  //           .slice(0, groupIndex)
  //           .reduce((x0, [_, sigs]) => x0 + sigs.length, 0) + i
  //     ),
  //     y: signatures.map((e) => e.contribution),
  //     //text: signatures.map((e, i) => e.mutationType),
  //     //hovertemplate: "%{signatures.map((e, i) => e.mutationType)}, %{y}",
  //     hoverinfo: 'x+y',
  //     showlegend: false,
  //   })
  // );
  const traces = data.map((group, groupIndex, array) => ({
    name: group.indel,
    type: 'bar',
    marker: { color: colors[group.indel].shape },
    x: [...group.data.keys()].map(
      (e) =>
        e +
        array
          .slice(0, groupIndex)
          .reduce((lastIndex, b) => lastIndex + b.data.length, 0)
    ),
    y: group.data.map((e) => e.mutations),
    customdata: group.data.map((e) => ({ mutationType: e.mutationType })),
    hovertemplate: '%{customdata.mutationType}, %{y}<extra></extra>',
    showlegend: false,
  }));

  const shapeAnnotations = data.map((group, groupIndex, array) => ({
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x:
      array
        .slice(0, groupIndex)
        .reduce((lastIndex, b) => lastIndex + b.data.length, 0) +
      (group.data.length - 1) * 0.5,
    y: 1.01,
    text: `<b>${group.indel.slice(-1)}</b>`,
    showarrow: false,
    font: {
      size: 14,
      color: colors[group.indel].text,
    },
    align: 'center',
  }));

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
    y: 0.88,
    text:
      '<b>' +
      sample +
      ': ' +
      totalMutations.toLocaleString(undefined) +
      ' indels</b>',
    showarrow: false,
    font: {
      size: 18,
    },
    align: 'center',
  };

  const topShapes = data.map((group, groupIndex, array) => ({
    type: 'rect',
    xref: 'x',
    yref: 'paper',
    x0: array
      .slice(0, groupIndex)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.4),
    x1: array
      .slice(0, groupIndex + 1)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.6),
    y0: 1.07,
    y1: 1.01,
    fillcolor: colors[group.indel].shape,
    line: {
      width: 0,
    },
  }));

  const bottomShapes = data.map((group, groupIndex, array) => ({
    type: 'rect',
    xref: 'x',
    yref: 'paper',
    x0: array
      .slice(0, groupIndex)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.4),
    x1: array
      .slice(0, groupIndex + 1)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.6),
    y0: -0.01,
    y1: -0.05,
    fillcolor: colors[group.indel].shape,
    line: {
      width: 0,
    },
  }));

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: 500,
    //width:1080,
    autosize: true,
    xaxis: {
      showticklabels: false,
      showline: true,
      tickfont: { size: 11 },
      tickmode: 'array',
      tickvals: indelNames.map((_, i) => i),
      ticktext: indelNames.map((e) => e.index),
      linecolor: '#E0E0E0',
      linewidth: 1,
      mirror: 'all',
    },
    yaxis: {
      title: 'Number of Idels',
      autorange: false,
      range: [0, maxMutation * 1.2],
      linecolor: 'black',
      linewidth: 1,
      mirror: true,
    },

    shapes: [...topShapes, ...bottomShapes],
    annotations: [
      ...shapeAnnotations,
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
