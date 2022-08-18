import { groupBy } from 'lodash';

export default function pcBetweenSamples(rawData, args) {
  console.log(rawData);
  console.log(args);
  const samples = args.sample.split(',');

  const groupBySample = groupBy(rawData, 'sample');
  const sample1 = groupBySample[samples[0]].flat();
  const sample2 = groupBySample[samples[1]].flat();
  console.log(sample1);
  console.log(sample2);

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

  const groupByMutation1 = sample1.reduce((acc, e, i) => {
    const mutationRegex = /^(.{7})/;
    const mutation = e.mutationType.match(mutationRegex)[1];
    acc[mutation] = acc[mutation] ? [...acc[mutation], e] : [e];
    return acc;
  }, {});

  console.log(groupByMutation1);
  const sample1data = Object.entries(groupByMutation1).map(
    ([mutation, data]) => ({
      mutation,
      data,
    })
  );
  console.log(sample1data);

  const groupByMutation2 = sample2.reduce((acc, e, i) => {
    const mutationRegex = /^(.{7})/;
    const mutation = e.mutationType.match(mutationRegex)[1];
    acc[mutation] = acc[mutation] ? [...acc[mutation], e] : [e];
    return acc;
  }, {});
  console.log(groupByMutation2);
  const sample2data = Object.entries(groupByMutation2).map(
    ([mutation, data]) => ({
      mutation,
      data,
    })
  );
  console.log(sample2data);
  // const indelNames = sample1data
  //   .map((indel) =>
  //     indel.mutation.substring(0, 5) === '2:Del' ||
  //     indel.mutation.substring(0, 5) === '3:Del' ||
  //     indel.mutation.substring(0, 5) === '4:Del'
  //       ? indel.mutation.substring(0, 1)
  //       : indel.mutation
  //   )
  //   .flat();
  const indelNames = sample1data
    .map((indel) =>
      indel.data.map((e) => ({
        indel: e.mutationType,
        index:
          e.mutationType.substring(2, 5) === 'Del'
            ? +e.mutationType.slice(-1) + 1
            : e.mutationType.slice(-1),
      }))
    )
    .flat();

  console.log(indelNames);
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
  console.log(xLabelAnnotation);
  console.log(sample1data);
  const trace1 = sample1data.map((group, groupIndex, array) => ({
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
    y: group.data.map((e) => e.mutations),
    groupdata: group.data,
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
    hoverinfo: 'x+y',
    yaxis: 'y3',
  }));
  const trace2 = sample2data.map((group, groupIndex, array) => ({
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
    y: group.data.map((e) => e.mutations),
    groupdata: group.data,
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
    hoverinfo: 'x+y',
    yaxis: 'y2',
  }));
  console.log(trace1);
  const traces = [...trace1, ...trace2];
  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: 500,
    grid: {
      rows: 3,
      column: 1,
    },
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
      autorange: true,
      linecolor: '#D3D3D3',
      linewidth: 1,
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      showgrid: true,
      gridcolor: '#F5F5F5',
      domain: [0, 0.33],
    },
    yaxis2: {
      autorange: true,
      //range: [0, maxMutations * 1.3],
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      domain: [0.34, 0.66],
    },
    yaxis3: {
      autorange: true,
      //range: [0, maxMutations * 1.3],
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      domain: [0.67, 1],
    },

    //shapes: [...topShapes, ...bottomShapes],
    annotations: [
      //   ...shapeAnnotations,
      ...xLabelAnnotation,
      //   ...annotationsIDTopLabel,
      //   ...annotationsIDBotLabel,
      //   sampleAnnotation,
    ],
  };

  return { traces, layout };
}
