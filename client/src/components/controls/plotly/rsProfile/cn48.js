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

  var sortOrder = ['0-100kb', '100kb-1Mb', '1Mb-10Mb', '10Mb-40Mb', '>40Mb']; // Declare a array that defines the order of the elements to be sorted.

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

  const groupbyfirst2 = rawData.reduce((acc, e, i) => {
    const names = e.mutationType.split(':');

    acc[names[0] + ':' + names[1]] = acc[names[0] + ':' + names[1]]
      ? [...acc[names[0] + ':' + names[1]], e]
      : [e];
    return acc;
  }, {});
  console.log(groupbyfirst2);
  const groupbyfirst2Data = Object.entries(groupbyfirst2).map(
    ([mutation, data]) => ({
      mutation,
      data,
    })
  );
  console.log(groupbyfirst2Data);
  const thesort = (arr) => {
    // first grab the obj that is not getting sorted
    let first = arr.shift();
    // add the sortable tag
    let narr = arr
      .map((a) => ({ ...a, tag: a.mutation.split(':')[1] }))
      .sort((a, b) => b.tag.localeCompare(a.tag)) // sort it
      .map((a) => {
        delete a.tag;
        return a;
      }); // remove the custom tag
    narr.unshift(first);
    return narr;
  };

  const sortGroupByFirst2Data = thesort(groupbyfirst2Data);

  // function compareMutations(a, b) {
  //   let [countA, nameA] = a.mutation.split(':');
  //   let [countB, nameB] = b.mutation.split(':');
  //   let nameComparison = nameA.localeCompare(nameB);
  //   let countComparison = parseInt(countA) - parseInt(countB);
  //   if (nameComparison == 0) {
  //     return countComparison;
  //   } else {
  //     return nameComparison;
  //   }
  // }
  // const sortArray = groupbyfirst2Data.sort(compareMutations);
  // const firstPart = sortArray.slice(0, 4);
  // const secondPart = sortArray.slice(4, 10);
  // const customSortArray = [...secondPart, ...firstPart];
  // console.log(customSortArray);
  // console.log(firstPart);
  // console.log(secondPart);
  const dataD = groupByClusterData
    .map((indel) => indel.data.map((e) => e))
    .flat();
  console.log(dataD);
  const mutationTypeNames = dataD.map((group, i) => ({
    mutationType: group.mutationType.split(':')[2],
    index: i,
  }));

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

  console.log(sortGroupByFirst2Data);
  const topShapes = sortGroupByFirst2Data.map((group, groupIndex, array) => ({
    group: group,
    name: group.mutation,
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
    fillcolor: colors[group.mutation.split(':')[0]],
    line: {
      width: 0,
    },
    showlegend: false,
  }));
  console.log(topShapes);

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: 500,
    autosize: true,
    xaxis: {
      showticklabels: true,
      showline: true,
      tickangle: -90,
      //tickfont: { size: 11 },
      tickmode: 'array',

      linecolor: 'black',
      linewidth: 1,
      mirror: 'all',
      tickvals: mutationTypeNames.map((_, i) => i),
      ticktext: mutationTypeNames.map((e) => e.mutationType),
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
    // annotations: [
    //   //...topShapeAnnitations,
    //   //topShapeClusterAnnotation,
    //   //topShapeNonClusterAnnotation,
    //   //sampleAnnotation,
    // ],
  };

  return { traces, layout };
}
