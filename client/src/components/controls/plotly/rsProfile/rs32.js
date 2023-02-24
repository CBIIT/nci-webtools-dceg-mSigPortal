import { groupBy } from 'lodash';
import { rs32Color } from '../../utils/colors';
export default function RS32(rawData, sample) {
  const colors = rs32Color;

  const totalMutations = rawData.reduce(
    (total, indel) => total + indel.contribution,
    0
  );
  const maxMutation = Math.max(...rawData.map((indel) => indel.contribution));

  var sortOrder = ['1-10Kb', '10-100Kb', '100Kb-1Mb', '1Mb-10Mb', '>10Mb']; // Declare a array that defines the order of the elements to be sorted.

  const clusterd = [];
  const nonClustered = [];
  rawData.map((e) => {
    if (e.mutationType.substring(0, 3) === 'non') {
      nonClustered.push(e);
    } else {
      clusterd.push(e);
    }
  });

  const groupByIndelCluster = clusterd.reduce((acc, e, i) => {
    const indel = e.mutationType.substring(0, 13);

    acc[indel] = acc[indel] ? [...acc[indel], e] : [e];
    return acc;
  }, {});

  const groupByIndelNonCluster = nonClustered.reduce((acc, e, i) => {
    const indel = e.mutationType.substring(0, 17);

    acc[indel] = acc[indel] ? [...acc[indel], e] : [e];
    return acc;
  }, {});
  const clusterGroup = Object.entries(groupByIndelCluster).map(
    ([indel, data]) => ({
      indel: indel,
      data: data,
    })
  );

  const nonClusterGroup = Object.entries(groupByIndelNonCluster).map(
    ([indel, data]) => ({
      indel: indel,
      data: data,
    })
  );

  const data = [...clusterGroup, ...nonClusterGroup];

  const sortGroup1 = clusterGroup.map((element, index, array) => ({
    mutation: element.mutation,
    data: element.data.sort(function (a, b) {
      return (
        sortOrder.indexOf(a.mutationType.split('_')[2]) -
        sortOrder.indexOf(b.mutationType.split('_')[2])
      );
    }),
  }));

  const sortedData1 = sortGroup1
    .map((indel) => indel.data.map((e) => e))
    .flat();

  const sortGroup2 = nonClusterGroup.map((element, index, array) => ({
    mutation: element.mutation,
    data: element.data.sort(function (a, b) {
      return (
        sortOrder.indexOf(a.mutationType.split('_')[2]) -
        sortOrder.indexOf(b.mutationType.split('_')[2])
      );
    }),
  }));

  const sortedData2 = sortGroup2
    .map((indel) => indel.data.map((e) => e))
    .flat();

  const sortData = [...sortedData1, ...sortedData2];
  const mutationTypeNames = sortData.map((group, i) => ({
    mutationType: !group.mutationType.split('_')[2]
      ? group.mutationType.split('_')[1]
      : group.mutationType.split('_')[2],
    index: i,
  }));

  const traces = sortData.map((group, groupIndex, array) => ({
    group: group,
    name: group.indel,
    type: 'bar',
    marker: {
      color: colors[group.mutationType],
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
  const separateLine = {
    type: 'line',
    xref: 'x',
    yref: 'paper',
    x0: 15.5,
    x1: 15.5,
    y0: 0,
    y1: 1,
    line: {
      color: '#808080',
      width: 1,
    },
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
      tickangle: -90,
      tickfont: { size: 11 },
      tickmode: 'array',

      linecolor: 'black',
      linewidth: 1,
      mirror: 'all',
      tickvals: mutationTypeNames.map((_, i) => i),
      ticktext: mutationTypeNames.map((e) => e.mutationType),
    },
    yaxis: {
      title: {
        text: '<b>Percentage(%)</b>',
        font: {
          family: 'Times New Roman',
          size: 18,
        },
      },
      autorange: false,
      range: [0, maxMutation * 1.25],
      tickformat: ',.1%',
      linecolor: 'black',
      linewidth: 1,
      mirror: true,
    },

    shapes: [...topShapes, topShapeCluster, topShapeNonluster, separateLine],
    annotations: [
      ...topShapeAnnitations,
      topShapeClusterAnnotation,
      topShapeNonClusterAnnotation,
      sampleAnnotation,
    ],
  };

  return { traces, layout };
}
