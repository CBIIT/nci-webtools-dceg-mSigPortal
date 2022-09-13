import { groupBy } from 'lodash';
export default function CN48(rawData, sample) {
  const colors = {
    '0:0-100kb': '#F0F8FF',
    '0:100kb-1Mb': '#787CE6',
    '0:>1Mb': '#0000CD',
    '1:0-100kb': '#EBEBEB',
    '1:100kb-1Mb': '#C5C5C5',
    '1:1Mb-10Mb': '#9F9F9F',
    '1:10Mb-40Mb': '#797979',
    '1:>40Mb': '#545454',
    '2:0-100kb': '#F5FFFA',
    '2:100kb-1Mb': '#C0E2C3',
    '2:1Mb-10Mb': '#8BC48E',
    '2:10Mb-40Mb': '#56A858',
    '2:>40Mb': '#228B22',
    '3-4:0-100kb': '#FFF0F5',
    '3-4:100kb-1Mb': '#DEBDEB',
    '3-4:1Mb-10Mb': '#BE8BE1',
    '3-4:10Mb-40Mb': '#9D58D7',
    '3-4:>40Mb': '#7D26CD',
    '5-8:0-100kb': '#FFFAF0',
    '5-8:100kb-1Mb': '#F2DCB3',
    '5-8:1Mb-10Mb': '#E6BF78',
    '5-8:10Mb-40Mb': '#D9A23C',
    '5-8:>40Mb': '#CD8500',
    '9+:0-100kb': '#FFE4E1',
    '9+:100kb-1Mb': '#E2ADBC',
    '9+:1Mb-10Mb': '#C47798',
    '9+:10Mb-40Mb': '#A84074',
    '9+:>40Mb': '#8B0A50',
    0: '#0000CD',
    1: '#545454',
    2: '#228B22',
    '3-4': '#7D26CD',
    '5-8': '#CD8500',
    '9+': '#8B0A50',
  };

  const totalMutations = rawData.reduce(
    (total, indel) => total + indel.contribution,
    0
  );
  const maxMutation = Math.max(...rawData.map((indel) => indel.contribution));
  const hd = [];
  const LOH = [];
  const het = [];
  rawData.map((e) => {
    const names = e.mutationType.split(':');
    if (names[1] === 'LOH') {
      LOH.push(e);
    } else if (names[1] === 'het') {
      het.push(e);
    } else {
      hd.push(e);
    }
  });
  console.log(hd);
  console.log(LOH);
  console.log(het);
  const groupByCluster = rawData.reduce((acc, e, i) => {
    const names = e.mutationType.split(':');

    acc[names[1]] = acc[names[1]] ? [...acc[names[1]], e] : [e];
    return acc;
  }, {});
  console.log(groupByCluster);

  const groupByClusterData = Object.entries(groupByCluster).map(
    ([mutation, data]) => ({
      mutation,
      data,
    })
  );
  console.log(groupByClusterData);

  const groupByIndelhd = hd.reduce((acc, e, i) => {
    const indel = e.mutationType.substring(0, 13);

    acc[indel] = acc[indel] ? [...acc[indel], e] : [e];
    return acc;
  }, {});
  console.log(groupByIndelhd);
  const groupByIndelLOH = LOH.reduce((acc, e, i) => {
    const indel = e.mutationType.substring(0, 17);

    acc[indel] = acc[indel] ? [...acc[indel], e] : [e];
    return acc;
  }, {});

  const groupByIndelHet = het.reduce((acc, e, i) => {
    const indel = e.mutationType.substring(0, 17);

    acc[indel] = acc[indel] ? [...acc[indel], e] : [e];
    return acc;
  }, {});
  console.log(groupByIndelLOH);
  const groupByIndelhdGroup = Object.entries(groupByIndelhd).map(
    ([indel, data]) => ({
      indel,
      data,
    })
  );

  const groupByIndelLOHGroup = Object.entries(groupByIndelLOH).map(
    ([indel, data]) => ({
      indel,
      data,
    })
  );
  const groupByIndelHetGroup = Object.entries(groupByIndelHet).map(
    ([indel, data]) => ({
      indel,
      data,
    })
  );

  const data = [
    ...groupByIndelhdGroup,
    ...groupByIndelLOHGroup,
    ...groupByIndelHetGroup,
  ];
  console.log(data);

  const dataD = data.map((indel) => indel.data.map((e) => e)).flat();
  console.log(dataD);
  const mutationTypeNames = dataD
    .map((group) =>
      group.data.map((e, i) => ({
        mutationType: e.mutationType.split(':')[2],
      }))
    )
    .flat();
  console.log(mutationTypeNames);

  const traces = dataD.map((group, groupIndex, array) => ({
    group: group,
    name: group.indel,
    color: group.mutationType.split(':')[0] + group.mutationType.split(':')[2],
    type: 'bar',
    marker: {
      color:
        colors[
          group.mutationType.split(':')[0] +
            ':' +
            group.mutationType.split(':')[2]
        ],
      line: {
        color: 'black',
        width: 1,
      },
    },
    x: [group.mutationType],
    y: [group.contribution],
    customdata: {
      mutationType: group.mutationType,
      contribution: group.contribution,
    },
    hoverinfo: 'x+y',
    showlegend: false,
  }));
  console.log(traces);
  const traces1 = data.map((group, groupIndex, array) => ({
    group: group,
    name: group.indel,
    type: 'bar',
    marker: { color: colors[group.mutationType] },
    x: [...group.data.keys()].map(
      (e) =>
        e +
        array
          .slice(0, groupIndex)
          .reduce((lastIndex, b) => lastIndex + b.data.length, 0)
    ),
    y: group.data.map((e) => e.contribution),
    showlegend: false,
  }));
  const topShapes = data.map((group, groupIndex, array) => ({
    group: group,
    name: group.indel,
    type: 'rect',
    xref: 'x',
    yref: 'paper',
    marker: {
      color:
        colors[
          group.indel.substring(group.indel.length - 3, group.indel.length)
        ],
    },
    x0: array
      .slice(0, groupIndex)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.4),
    x1: array
      .slice(0, groupIndex + 1)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.6),
    y0: 1.07,
    y1: 1.01,
    fillcolor:
      colors[group.indel.substring(group.indel.length - 3, group.indel.length)],
    line: {
      width: 0,
    },
    showlegend: false,
  }));
  console.log(topShapes);

  const topShapeAnnitations = data.map((group, groupIndex, array) => ({
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
    text:
      group.indel.substring(group.indel.length - 3, group.indel.length) !==
      'tra'
        ? group.indel
            .substring(group.indel.length - 3, group.indel.length)
            .charAt(0)
            .toUpperCase() +
          group.indel
            .substring(group.indel.length - 3, group.indel.length)
            .slice(1)
        : 'T',
    showarrow: false,
    font: {
      size: 14,
      color: 'white',
    },
    align: 'center',
  }));
  const topShapeCluster = {
    type: 'rect',
    xref: 'x',
    yref: 'paper',
    x0: -0.4,
    x1: 15.4,
    y0: 1.14,
    y1: 1.08,
    fillcolor: '#808080',
    line: {
      width: 0,
    },
  };
  const topShapeClusterAnnotation = {
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: 7.5,
    y: 1.08,
    text: `Clustered`,
    showarrow: false,
    font: {
      size: 14,
      color: 'white',
    },
    align: 'center',
  };
  const topShapeNonluster = {
    type: 'rect',
    xref: 'x',
    yref: 'paper',
    x0: 15.6,
    x1: 31.4,
    y0: 1.14,
    y1: 1.08,
    fillcolor: '#000000',
    line: {
      width: 0,
    },
  };
  const topShapeNonClusterAnnotation = {
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: 23.5,
    y: 1.08,
    text: `Non-Clustered`,
    showarrow: false,
    font: {
      size: 14,
      color: 'white',
    },
    align: 'center',
  };

  const sampleAnnotation = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: 0.01,
    y: 0.88,
    text: '<b>' + sample + '</b>',
    showarrow: false,
    font: {
      size: 18,
      family: 'Arial',
    },
    align: 'center',
  };
  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: 500,
    autosize: true,
    xaxis: {
      showticklabels: true,
      showline: true,
      //tickangle: -90,
      //tickfont: { size: 11 },
      tickmode: 'array',

      linecolor: 'black',
      linewidth: 1,
      mirror: 'all',
      // tickvals: mutationTypeNames.map((_, i) => i),
      // ticktext: mutationTypeNames.map((e) => e.mutationType),
    },
    yaxis: {
      title: {
        text: '<b>Proportion</b>',
        font: {
          family: 'Times New Roman',
          size: 18,
        },
      },
      autorange: false,
      range: [0, maxMutation * 1.25],
      //tickformat: ',.1%',
      linecolor: 'black',
      linewidth: 1,
      mirror: true,
    },

    shapes: [...topShapes],
    annotations: [
      ...topShapeAnnitations,
      topShapeClusterAnnotation,
      topShapeNonClusterAnnotation,
      sampleAnnotation,
    ],
  };

  return { traces, layout };
}
