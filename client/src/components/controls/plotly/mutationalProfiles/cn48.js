import { getMaxMutations, createSampleAnnotation } from './utils.js';

export default function CN48(apiData, title = '') {
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
    homdel: '#0000CD',
    LOH: '#FFFFFF',
    het: 'black',
  };

  const maxMutation = getMaxMutations(apiData);

  var sortOrder = [
    '0-100kb',
    '100kb-1Mb',
    '>1Mb',
    '1Mb-10Mb',
    '10Mb-40Mb',
    '>40Mb',
  ]; // Declare a array that defines the order of the elements to be sorted.

  const groupByCluster = apiData.reduce((acc, e, i) => {
    const names = e.mutationType.split(':');

    acc[names[1]] = acc[names[1]] ? [...acc[names[1]], e] : [e];
    return acc;
  }, {});

  const groupByClusterData = Object.entries(groupByCluster).map(
    ([mutation, data]) => ({
      mutation: mutation,
      data: data,
    })
  );

  const groupbyfirst2 = apiData.reduce((acc, e, i) => {
    const names = e.mutationType.split(':');

    acc[names[0] + ':' + names[1]] = acc[names[0] + ':' + names[1]]
      ? [...acc[names[0] + ':' + names[1]], e]
      : [e];
    return acc;
  }, {});

  const groupbyfirst2Data = Object.entries(groupbyfirst2).map(
    ([mutation, data]) => ({
      mutation,
      data,
    })
  );

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

  const sortGroupByFirst2DataInside = sortGroupByFirst2Data.map(
    (element, index, array) => ({
      mutation: element.mutation,
      data: element.data.sort(function (a, b) {
        return (
          sortOrder.indexOf(a.mutationType.split(':')[2]) -
          sortOrder.indexOf(b.mutationType.split(':')[2])
        );
      }),
    })
  );

  const sortedData = sortGroupByFirst2DataInside
    .map((group) => group.data.map((e) => e))
    .flat();
  const mutationTypeNames = sortedData.map((group, i) => ({
    mutationType: group.mutationType.split(':')[2],
    index: i,
  }));

  const traces = sortedData.map((group, groupIndex, array) => ({
    group: group,
    name: group.mutationType,
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
    customdata: [
      {
        type: group.mutationType.split(':')[1],
        xval: group.mutationType.split(':')[2],
        contribution: group.contribution,
      },
    ],
    //hoverinfo: 'x+y',
    hovertemplate:
      '<b>%{customdata.type}</b><br>' +
      '%{customdata.xval} <br>Proportion: %{y}<extra></extra>',
    showlegend: false,
  }));

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
      width: 1,
    },
    showlegend: false,
  }));
  const topShapeAnnitations = sortGroupByFirst2Data.map(
    (group, groupIndex, array) => ({
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
      text: group.mutation.split(':')[0],
      showarrow: false,
      font: {
        size: 14,
        color: 'white',
      },
      align: 'center',
    })
  );
  const topTitleShapes = groupByClusterData.map((group, groupIndex, array) => ({
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
    y0: 1.14,
    y1: 1.08,
    fillcolor: colors[group.mutation.split(':')[0]],
    line: {
      width: 1,
    },
    showlegend: false,
  }));

  const topTitleShapesAnnitations = groupByClusterData.map(
    (group, groupIndex, array) => ({
      group: group,
      xref: 'x',
      yref: 'paper',
      xanchor: 'bottom',
      yanchor: 'bottom',
      x:
        array
          .slice(0, groupIndex)
          .reduce((lastIndex, b) => lastIndex + b.data.length, 0) +
        (group.data.length - 1) * 0.5,
      y: 1.08,
      text: group.mutation.charAt(0).toUpperCase() + group.mutation.slice(1),
      showarrow: false,
      font: {
        size: 14,
        color: group.mutation === 'LOH' ? 'black' : 'white',
      },
      align: 'center',
    })
  );
  const sampleAnnotation = createSampleAnnotation(apiData);
  const layout = {
    title: `<b>${title}</b>`,
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

    shapes: [...topShapes, ...topTitleShapes],
    annotations: [
      ...topShapeAnnitations,
      ...topTitleShapesAnnitations,
      sampleAnnotation,
    ],
  };

  return { traces, layout };
}
