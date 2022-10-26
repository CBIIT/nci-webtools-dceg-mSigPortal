import {
  getMaxMutations,
  getTotalMutations,
  groupDataByMutation,
} from './utils';

export default function ID83(apiData) {
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

  const indelRegex = /^(.{7})/;
  const mutationGroupSort = (a, b) => {
    const order = Object.keys(colors);
    return order.indexOf(a.mutation) - order.indexOf(b.mutation);
  };

  const data = groupDataByMutation(apiData, indelRegex, mutationGroupSort);

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
      'Microhomology Length',
    ],
    arrayIDAnnXLabel = [5, 18.5, 35, 60, 76];

  const totalMutations = getTotalMutations(apiData);
  const maxMutation = getMaxMutations(apiData);

  const indelNames = data
    .map((group) =>
      group.data.map((e) => ({
        indel: group.mutation,
        index:
          group.mutation.includes('Del') && group.mutation.slice(-1) != 'M'
            ? +e.mutationType.slice(-1) + 1
            : e.mutationType.slice(-1),
      }))
    )
    .flat();

  const traces = data.map((group, groupIndex, array) => ({
    name: group.mutation,
    type: 'bar',
    marker: { color: colors[group.mutation].shape },
    x: [...group.data.keys()].map(
      (e) =>
        e +
        array
          .slice(0, groupIndex)
          .reduce((lastIndex, b) => lastIndex + b.data.length, 0)
    ),
    y: group.data.map((e) => e.mutations || e.contribution),
    groupdata: group.data,
    //customdata: group.data.map((e) => ({ mutationType: e.mutationType })),
    customdata: group.data.map((e) => ({
      mutationOrder: e.mutationType.substring(0, 1),
      mutationType:
        e.mutationType.substring(2, 5) === 'Del' ? 'Deletion' : 'Insertion',
      extraValue: e.mutationType.substring(6, 7),
      xval:
        e.mutationType.substring(2, 5) === 'Del'
          ? +e.mutationType.slice(-1) + 1
          : e.mutationType.slice(-1),
    })),
    hovertemplate:
      '<b>%{customdata.mutationOrder} bp %{customdata.mutationType}, %{customdata.extraValue}, %{customdata.xval}</b><br>' +
      '%{y} indels<extra></extra>',
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
    text: `<b>${
      group.mutation[0] == '1' ? group.mutation.slice(-1) : group.mutation[0]
    }</b>`,
    showarrow: false,
    font: {
      size: 14,
      color: colors[group.mutation].text,
    },
    align: 'center',
  }));

  const xLabelAnnotation = indelNames.map((indel, index) => ({
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: index,
    y: -0.1,
    text: '<b>' + indel.index + '</b>',
    showarrow: false,
    font: {
      size: 12,
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
      size: 16,
      family: 'Times New Roman',
    },
    align: 'center',
  }));

  const sampleAnnotation = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: 0.01,
    y: 0.88,
    text: apiData[0].sample
      ? `<b>${apiData[0].sample}: ${totalMutations.toLocaleString()} Indels</b>`
      : `<b>${apiData[0].signatureName}</b>`,
    showarrow: false,
    font: {
      size: 24,
      family: 'Arial',
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
    fillcolor: colors[group.mutation].shape,
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
    fillcolor: colors[group.mutation].shape,
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
      linecolor: 'black',
      linewidth: 1,
      mirror: 'all',
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
      range: [0, maxMutation * 1.2],
      linecolor: 'black',
      linewidth: 1,
      mirror: true,
    },

    shapes: [...topShapes, ...bottomShapes],
    annotations: [
      ...shapeAnnotations,
      ...xLabelAnnotation,
      ...annotationsIDTopLabel,
      ...annotationsIDBotLabel,
      sampleAnnotation,
    ],
  };

  return { traces, layout };
}
