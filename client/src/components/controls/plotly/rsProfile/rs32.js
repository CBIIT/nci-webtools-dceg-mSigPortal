import { groupBy } from 'lodash';
export default function RS32(rawData, sample) {
  const colors = {
    'clustered_del_>10Mb': 'deeppink',
    'non-clustered_del_>10Mb': 'deeppink',
    'clustered_del_1Mb-10Mb': 'hotpink',
    'non-clustered_del_1Mb-10Mb': 'hotpink',
    'clustered_del_10-100Kb': 'lightpink',
    'non-clustered_del_10-100Kb': 'lightpink',
    'clustered_del_100Kb-1Mb': 'palevioletred',
    'non-clustered_del_100Kb-1Mb': 'palevioletred',
    'clustered_del_1-10Kb': 'lavenderblush',
    'non-clustered_del_1-10Kb': 'lavenderblush',
    'clustered_tds_>10Mb': 'saddlebrown',
    'non-clustered_tds_>10Mb': 'saddlebrown',
    'clustered_tds_1Mb-10Mb': 'sienna',
    'non-clustered_tds_1Mb-10Mb': 'sienna',
    'clustered_tds_10-100Kb': 'sandybrown',
    'non-clustered_tds_10-100Kb': 'sandybrown',
    'clustered_tds_100Kb-1Mb': 'peru',
    'non-clustered_tds_100Kb-1Mb': 'peru',
    'clustered_tds_1-10Kb': 'linen',
    'non-clustered_tds_1-10Kb': 'linen',
    'clustered_inv_>10Mb': 'rebeccapurple',
    'non-clustered_inv_>10Mb': 'rebeccapurple',
    'clustered_inv_1Mb-10Mb': 'blueviolet',
    'non-clustered_inv_1Mb-10Mb': 'blueviolet',
    'clustered_inv_10-100Kb': 'plum',
    'non-clustered_inv_10-100Kb': 'plum',
    'clustered_inv_100Kb-1Mb': 'mediumorchid',
    'non-clustered_inv_100Kb-1Mb': 'mediumorchid',
    'clustered_inv_1-10Kb': 'thistle',
    'non-clustered_inv_1-10Kb': 'thistle',
    'clustered_trans': 'thistle',
    'non-clustered_trans': 'thistle',
  };

  const clusterd = [];
  const nonClustered = [];
  rawData.map((e) => {
    if (e.mutationType.substring(0, 3) === 'non') {
      nonClustered.push(e);
    } else {
      clusterd.push(e);
    }
  });

  console.log(nonClustered);
  console.log(clusterd);

  const groupByIndelCluster = clusterd.reduce((acc, e, i) => {
    const indel = e.mutationType.substring(0, 13);

    acc[indel] = acc[indel] ? [...acc[indel], e] : [e];
    return acc;
  }, {});
  console.log(groupByIndelCluster);

  const indelNames = clusterd
    .map((e, i) => ({
      indel: e.mutationType.substring(14, e.mutationType.length - 1),
      index: i,
    }))
    .flat();
  console.log(indelNames);

  const arrayTop1 = ['Clustered', 'Non-Clustered'];
  const arrayTop2 = [''];

  const traces = rawData.map((group, groupIndex, array) => ({
    group: group,
    name: 'Clustered',
    type: 'bar',
    marker: {
      color: colors[group.mutationType],
      line: {
        color: 'black',
        width: 1.5,
      },
    },
    x: [group.mutationType],
    y: [group.contribution],
    hoverinfo: 'x+y',
    showlegend: false,
  }));
  console.log(traces);

  const topShapes = rawData.map((group, groupIndex, array) => ({
    type: 'rect',
    xref: 'x',
    yref: 'paper',
    x0: array
      .slice(0, groupIndex)
      .reduce((lastIndex, e) => lastIndex + e.length, -0.4),
    x1: array
      .slice(0, groupIndex + 1)
      .reduce((lastIndex, e) => lastIndex + e.length, -0.6),
    y0: 1.07,
    y1: 1.01,
    fillcolor: colors[group.mutationType],
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
      showticklabels: true,
      showline: true,
      tickfont: { size: 11 },
      tickmode: 'array',
      //tickvals: indelNames.map((_, i) => i),
      //ticktext: indelNames.map((e) => e.index),
      linecolor: 'black',
      linewidth: 1,
      mirror: 'all',
    },
    yaxis: {
      title: {
        text: '<b>Percentage(%)</b>',
        font: {
          family: 'Times New Roman',
          size: 18,
        },
      },
      autorange: true,
      //range: [0, maxMutation * 1.2],
      linecolor: 'black',
      linewidth: 1,
      mirror: true,
    },

    shapes: [...topShapes],
    // annotations: [
    //   ...shapeAnnotations,
    //   ...xLabelAnnotation,
    //   ...annotationsIDTopLabel,
    //   ...annotationsIDBotLabel,
    //   sampleAnnotation,
    // ],
  };

  return { traces, layout };
}
