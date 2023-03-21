import { createSampleAnnotation } from './utils.js';

export default function ID415(data, title = '') {
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

  const transcribed = data.filter((e) => /^T:/.test(e.mutationType));
  const untranscribed = data.filter((e) => /^U:/.test(e.mutationType));

  const totalMutations =
    transcribed.reduce((total, e) => total + e.mutations, 0) +
    untranscribed.reduce((total, e) => total + e.mutations, 0);

  const maxMutation = Math.max(
    ...[
      ...transcribed.map((e) => e.mutations),
      ...untranscribed.map((e) => e.mutations),
    ]
  );

  ///// --------- T Group ------------///////
  const T_groupByMutation = transcribed.reduce((groups, e, i) => {
    const mutationRegex = /^.{2,9}/;
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

  const T_groupByFirstGroup = Object.fromEntries(
    Object.entries(T_groupByMutation).slice(0, 4)
  );

  const T_groupByMutationID = transcribed.reduce((groups, e) => {
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

  const T_groupR = T_groupByMutationID['R'].reduce((r, a) => {
    let m;
    m = a.mutationType.match(a.mutationType.substr(4, 3));
    const s = {
      mutationType: a.mutationType,
      contribution: a.contribution,
    };
    r[m] = r[m] ? [...r[m], a] : [s];
    return r;
  }, {});

  const T_groupRDel = T_groupR['Del'].reduce((r, a) => {
    let m;
    m = a.mutationType.match(a.mutationType.substr(0, 7));
    const s = {
      mutationType: a.mutationType,
      contribution: a.contribution,
    };
    r[m] = r[m] ? [...r[m], a] : [s];
    return r;
  }, {});

  const T_groupRIns = T_groupR['Ins'].reduce((r, a) => {
    let m;
    m = a.mutationType.match(a.mutationType.substr(0, 7));
    const s = {
      mutationType: a.mutationType,
      contribution: a.contribution,
    };
    r[m] = r[m] ? [...r[m], a] : [s];
    return r;
  }, {});

  const T_groupM = T_groupByMutationID['M'].reduce((r, a) => {
    let m;
    m = a.mutationType.match(a.mutationType.substr(0, 7));
    const s = {
      mutationType: a.mutationType,
      contribution: a.contribution,
    };
    r[m] = r[m] ? [...r[m], a] : [s];
    return r;
  }, {});
  const T_arrayID1 = Object.keys(T_groupByFirstGroup).map(function (key) {
    return T_groupByFirstGroup[key];
  });
  const T_arrayID2 = Object.keys(T_groupRDel).map(function (key) {
    return T_groupRDel[key];
  });
  const T_arrayID3 = Object.keys(T_groupRIns).map(function (key) {
    return T_groupRIns[key];
  });
  const T_arrayID4 = Object.keys(T_groupM).map(function (key) {
    return T_groupM[key];
  });

  const T_arrayID = [
    ...T_arrayID1,
    ...T_arrayID2,
    ...T_arrayID3,
    ...T_arrayID4,
  ];

  const T_flatSorted = Object.values(T_arrayID).flat();

  ///// --------- U Group ------------///////
  const U_groupByMutation = untranscribed.reduce((groups, e, i) => {
    const mutationRegex = /^.{2,9}/;
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

  const U_groupByFirstGroup = Object.fromEntries(
    Object.entries(U_groupByMutation).slice(0, 4)
  );

  const U_groupByMutationID = untranscribed.reduce((groups, e) => {
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

  const U_groupR = U_groupByMutationID['R'].reduce((r, a) => {
    let m;
    m = a.mutationType.match(a.mutationType.substr(4, 3));
    const s = {
      mutationType: a.mutationType,
      contribution: a.contribution,
    };
    r[m] = r[m] ? [...r[m], a] : [s];
    return r;
  }, {});

  const U_groupRDel = U_groupR['Del'].reduce((r, a) => {
    let m;
    m = a.mutationType.match(a.mutationType.substr(0, 7));
    const s = {
      mutationType: a.mutationType,
      contribution: a.contribution,
    };
    r[m] = r[m] ? [...r[m], a] : [s];
    return r;
  }, {});

  const U_groupRIns = U_groupR['Ins'].reduce((r, a) => {
    let m;
    m = a.mutationType.match(a.mutationType.substr(0, 7));
    const s = {
      mutationType: a.mutationType,
      contribution: a.contribution,
    };
    r[m] = r[m] ? [...r[m], a] : [s];
    return r;
  }, {});

  const U_groupM = U_groupByMutationID['M'].reduce((r, a) => {
    let m;
    m = a.mutationType.match(a.mutationType.substr(0, 7));
    const s = {
      mutationType: a.mutationType,
      contribution: a.contribution,
    };
    r[m] = r[m] ? [...r[m], a] : [s];
    return r;
  }, {});

  const U_arrayID1 = Object.keys(U_groupByFirstGroup).map(function (key) {
    return U_groupByFirstGroup[key];
  });
  const U_arrayID2 = Object.keys(U_groupRDel).map(function (key) {
    return U_groupRDel[key];
  });
  const U_arrayID3 = Object.keys(U_groupRIns).map(function (key) {
    return U_groupRIns[key];
  });
  const U_arrayID4 = Object.keys(U_groupM).map(function (key) {
    return U_groupM[key];
  });

  const U_arrayID = [
    ...U_arrayID1,
    ...U_arrayID2,
    ...U_arrayID3,
    ...U_arrayID4,
  ];
  const U_flatSorted = Object.values(U_arrayID).flat();

  //// ----------- plot ------------------//
  const tracesT = {
    name: 'Transcrribed Strand',
    type: 'bar',
    marker: { color: '#004765' },
    customedata: T_flatSorted.map((e, i, a) => ({
      mutationOrder: e.mutationType.substring(0, 1),
      mutationType:
        e.mutationType.substring(2, 5) === 'Del' ? 'Deletion' : 'Insertion',
      extraValue: e.mutationType.substring(6, 7),
      xval:
        e.mutationType.substring(2, 5) === 'Del'
          ? +e.mutationType.slice(-1) + 1
          : e.mutationType.slice(-1),
    })),
    x: T_flatSorted.map((element, index, array) => index),
    y: T_flatSorted.map((element, index, array) => element.contribution),
    hovertemplate: '<b>Transcrribed Strand</b><br>%{y} indels <extra></extra>',
    //hoverinfo: 'x+y',
    showlegend: true,
  };

  const tracesU = {
    name: 'Untranscribed Strand',
    type: 'bar',
    marker: { color: '#E32925' },
    x: U_flatSorted.map((element, index, array) => index),
    y: U_flatSorted.map((element, index, array) => element.contribution),
    hovertemplate: '<b>Untranscribed Strand</b><br>%{y} indels <extra></extra>',
    //hoverinfo: 'x+y',
    showlegend: true,
  };

  Object.values(T_arrayID).forEach((group) => {
    if (group.length > 1) {
      arrayIDAnnotationTop.push(
        group[Math.floor(group.length / 2)].mutationType
      );
    } else {
      arrayIDAnnotationTop.push(group[0].mutationType);
    }
    group.forEach((e) => {
      let lastNum = e.mutationType.substring(
        e.mutationType.length - 1,
        e.mutationType.length
      );
      let newNum;
      if (
        e.mutationType.substring(4, 9) === 'Del:C' ||
        e.mutationType.substring(4, 9) === 'Del:T' ||
        e.mutationType.substring(4, 9) === 'Del:R'
      ) {
        lastNum = +lastNum + 1;
      }
      if (
        (e.mutationType.substring(4, 9) === 'Del:C' ||
          e.mutationType.substring(4, 9) === 'Del:T' ||
          e.mutationType.substring(4, 9) === 'Del:R') &
        (+lastNum > 5)
      ) {
        newNum = lastNum + '+';
      } else if (
        e.mutationType.substring(4, 9) !== 'Del:C' &&
        e.mutationType.substring(4, 9) !== 'Del:T' &&
        e.mutationType.substring(4, 9) !== 'Del:R' &&
        +lastNum > 4
      ) {
        newNum = lastNum + '+';
      } else {
        newNum = lastNum;
      }
      arrayIDAnnotationBot.push(newNum);
    });
  });

  const traces = [tracesT, tracesU];

  const annotations1 = Object.entries(T_arrayID).map(
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
          : `<b>${signatures[0].mutationType.substring(2, 3)}</b>`,
      showarrow: false,
      font: {
        size: 14,
        color:
          annotationColors[
            signatures[0].mutationType.substring(
              2,
              signatures[0].mutationType.length - 2
            )
          ],
      },
      align: 'center',
    })
  );

  const annotations2 = arrayIDAnnotationBot.map((num, index) => ({
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: index,
    y: -0.1,
    text: '<b>' + num + '</b>',
    showarrow: false,
    font: {
      size: 12,
      family: 'Times New Roman',
    },
    align: 'center',
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
      size: 16,
      family: 'Times New Roman',
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
      size: 15,
      family: 'Times New Roman',
    },
    align: 'center',
  }));

  const sampleAnnotation = createSampleAnnotation(data);

  const shapes1 = Object.entries(T_arrayID).map(
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
            2,
            signatures[0].mutationType.length - 2
          )
        ],
      line: {
        width: 0,
      },
    })
  );

  const shapes2 = Object.entries(T_arrayID).map(
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
            2,
            signatures[0].mutationType.length - 2
          )
        ],
      line: {
        width: 0,
      },
    })
  );

  const layout = {
    title: `<b>${title}</b>`,
    hoverlabel: { bgcolor: '#FFF' },
    height: 500,
    width: 1100,
    autosize: false,
    legend: {
      x: 1,
      xanchor: 'right',
      y: 1,
      bgcolor: '#FFFFFF',
      bordercolor: '#E1E1E1',
      borderwidth: 1,
    },
    xaxis: {
      showticklabels: false,
      showline: true,
      tickangle: -90,
      tickfont: {
        size: 10,
      },
      tickmode: 'array',
      tickvals: T_flatSorted.map((_, i) => i),
      ticktext: T_flatSorted.map((e) => e.mutationType),
      linecolor: 'black',
      linewidth: 1,
      mirror: true,
    },
    yaxis: {
      title: {
        text: '<b>Number of Indels</b>',
        font: {
          family: 'Times New Roman',
          size: 18,
        },
      },
      autorange: false,
      range: [0, maxMutation + maxMutation * 0.2],
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

  return { traces, layout };
}
