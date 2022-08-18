import { faExternalLinkAlt } from '@fortawesome/free-solid-svg-icons';
import { groupBy } from 'lodash';

export default function pcBetweenSamples(rawData, args) {
  const samples = args.sample.split(',');
  const groupBySample = groupBy(rawData, 'sample');
  const sample1 = groupBySample[samples[0]].flat();
  const sample2 = groupBySample[samples[1]].flat();
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
  const mutationRegex = /^(.{7})/;
  const groupByMutation1 = sample1.reduce((acc, e, i) => {
    const mutation = e.mutationType.match(mutationRegex)[1];
    acc[mutation] = acc[mutation] ? [...acc[mutation], e] : [e];
    return acc;
  }, {});

  const sample1data = Object.entries(groupByMutation1).map(
    ([mutation, data]) => ({
      mutation,
      data,
    })
  );

  const groupByMutation2 = sample2.reduce((acc, e, i) => {
    const mutation = e.mutationType.match(mutationRegex)[1];
    acc[mutation] = acc[mutation] ? [...acc[mutation], e] : [e];
    return acc;
  }, {});
  const sample2data = Object.entries(groupByMutation2).map(
    ([mutation, data]) => ({
      mutation,
      data,
    })
  );

  const totalMutations1 = sample1data.reduce(
    (total, mutation) =>
      total +
      mutation.data.reduce((mutationSum, e) => mutationSum + e.mutations, 0),
    0
  );
  const maxMutation1 = Math.max(
    ...sample1data
      .map((mutation) =>
        mutation.data.map((e) => e.mutations / totalMutations1)
      )
      .flat()
  );

  const totalMutations2 = sample2data.reduce(
    (total, mutation) =>
      total +
      mutation.data.reduce((mutationSum, e) => mutationSum + e.mutations, 0),
    0
  );
  const maxMutation2 = Math.max(
    ...sample2data
      .map((mutation) =>
        mutation.data.map((e) => e.mutations / totalMutations2)
      )
      .flat()
  );
  const maxMutations = Math.max(maxMutation1, maxMutation2);
  const group1 = groupBy(sample1, 'mutationType');
  Object.keys(group1);
  const group2 = groupBy(sample2, 'mutationType');
  Object.keys(group2);

  let sampleDifferences = [];
  let s1mutations = [];
  let s2mutations = [];

  for (let mutationType of Object.keys(group1)) {
    const a = group1[mutationType][0];
    const b = group2[mutationType][0];
    const mutations =
      a.mutations / totalMutations1 - b.mutations / totalMutations2;
    sampleDifferences.push({ mutationType, mutations });
    s1mutations.push(a.mutations / totalMutations1);
    s2mutations.push(b.mutations / totalMutations2);
  }
  const groupByMutation3 = groupBy(
    sampleDifferences,
    (s) => s.mutationType.match(mutationRegex)[1]
  );
  const sample3data = Object.entries(groupByMutation3).map(
    ([mutation, data]) => ({
      mutation,
      data,
    })
  );

  const squarediff = sampleDifferences.map((e) => Math.pow(e.mutations, 2));
  const rss = squarediff.reduce((a, b, i) => a + b, 0).toExponential(3);

  function dotp(x, y) {
    function dotp_sum(a, b) {
      return a + b;
    }
    function dotp_times(a, i) {
      return x[i] * y[i];
    }
    return x.map(dotp_times).reduce(dotp_sum, 0);
  }

  function cosineSimilarity(A, B) {
    var similarity =
      dotp(A, B) / (Math.sqrt(dotp(A, A)) * Math.sqrt(dotp(B, B)));
    return similarity;
  }
  const cosine = cosineSimilarity(s1mutations, s2mutations).toFixed(3);

  const indelNames0 = sample1data.map((indel) =>
    indel.data.map((e) => ({
      indel: e.mutationType,
      index:
        e.mutationType.substring(2, 5) === 'Del' &&
        e.mutationType.substring(2, 7) !== 'Del:M'
          ? +e.mutationType.slice(-1) + 1
          : e.mutationType.slice(-1),
    }))
  );
  const indelNames = indelNames0
    .map((e) =>
      e.map((el, i) => ({
        indel: el.indel,
        index:
          i < e.length - 1
            ? el.index
            : e.length >= 5
            ? el.index + '+'
            : el.index,
        length: e.length,
        i: i,
      }))
    )
    .flat();

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
    y: group.data.map((e) => e.mutations / totalMutations1),
    customdata: group.data.map((e) => ({
      mutationType: e.mutationType.substring(0, 7),
      xval:
        e.mutationType.substring(2, 5) === 'Del' &&
        e.mutationType.substring(2, 7) !== 'Del:M'
          ? +e.mutationType.slice(-1) + 1
          : e.mutationType.slice(-1),
    })),

    hovertemplate:
      '<b>%{customdata.mutationType}, %{customdata.xval}</b><br>' +
      '%{y}<extra></extra>',
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
    y: group.data.map((e) => e.mutations / totalMutations2),
    customdata: group.data.map((e) => ({
      mutationType: e.mutationType.substring(0, 7),
      xval:
        e.mutationType.substring(2, 5) === 'Del' &&
        e.mutationType.substring(2, 7) !== 'Del:M'
          ? +e.mutationType.slice(-1) + 1
          : e.mutationType.slice(-1),
    })),

    hovertemplate:
      '<b>%{customdata.mutationType}, %{customdata.xval}</b><br>' +
      '%{y}<extra></extra>',
    showlegend: false,
    hoverinfo: 'x+y',
    yaxis: 'y2',
  }));
  const trace3 = sample3data.map((group, groupIndex, array) => ({
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
    customdata: group.data.map((e) => ({
      mutationType: e.mutationType.substring(0, 7),
      xval:
        e.mutationType.substring(2, 5) === 'Del' &&
        e.mutationType.substring(2, 7) !== 'Del:M'
          ? +e.mutationType.slice(-1) + 1
          : e.mutationType.slice(-1),
    })),

    hovertemplate:
      '<b>%{customdata.mutationType}, %{customdata.xval}</b><br>' +
      '%{y}<extra></extra>',
    showlegend: false,
    hoverinfo: 'x+y',
  }));
  const traces = [...trace1, ...trace2, ...trace3];

  const shapeTop = sample1data.map((group, groupIndex, array) => ({
    type: 'rect',
    xref: 'x',
    yref: 'paper',
    x0: array
      .slice(0, groupIndex)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
    x1: array
      .slice(0, groupIndex + 1)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
    y0: 1.0,
    y1: 1.05,
    fillcolor: colors[group.mutation].shape,
    line: {
      width: 1,
    },
  }));

  const shapeLine3 = sample1data.map((group, groupIndex, array) => ({
    type: 'rect',
    xref: 'x',
    yref: 'paper',
    x0: array
      .slice(0, groupIndex)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
    x1: array
      .slice(0, groupIndex + 1)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
    y0: 0.67,
    y1: 1,

    line: {
      width: 1,
    },
  }));

  const shapeLine2 = sample2data.map((group, groupIndex, array) => ({
    type: 'rect',
    xref: 'x',
    yref: 'paper',
    x0: array
      .slice(0, groupIndex)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
    x1: array
      .slice(0, groupIndex + 1)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
    y0: 0.34,
    y1: 0.66,

    line: {
      width: 1,
    },
  }));

  const shapeLine1 = sample3data.map((group, groupIndex, array) => ({
    type: 'rect',
    xref: 'x',
    yref: 'paper',
    x0: array
      .slice(0, groupIndex)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
    x1: array
      .slice(0, groupIndex + 1)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
    y0: 0.33,
    y1: 0,

    line: {
      width: 1,
    },
  }));
  const shapeRight3 = {
    type: 'rect',
    xref: 'paper',
    yref: 'paper',
    x0: 1,
    x1: 1.02,
    y0: 0.67,
    y1: 1,
    fillcolor: '#F0F0F0',
  };
  const shapeRight2 = {
    type: 'rect',
    xref: 'paper',
    yref: 'paper',
    x0: 1,
    x1: 1.02,
    y0: 0.34,
    y1: 0.66,
    fillcolor: '#F0F0F0',
  };

  const shapeRight1 = {
    type: 'rect',
    xref: 'paper',
    yref: 'paper',
    x0: 1,
    x1: 1.02,
    y0: 0,
    y1: 0.33,
    fillcolor: '#F0F0F0',
  };

  const annotationLabelRight3 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'middle',
    yanchor: 'middle',
    align: 'center',
    x: 1.0175,
    y: 0.835,
    text: samples[0],
    textangle: 90,
    showarrow: false,
  };

  const annotationLabelRight2 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'middle',
    yanchor: 'middle',
    align: 'center',
    x: 1.0175,
    y: 0.505,
    text: samples[1],
    textangle: 90,
    showarrow: false,
  };

  const annotationLabelRight1 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'middle',
    yanchor: 'middle',
    align: 'center',
    x: 1.0175,
    y: 0.165,
    text: 'Difference',
    textangle: 90,
    showarrow: false,
  };
  const mutationAnnotation = sample1data.map((group, groupIndex, array) => ({
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
      group.mutation.substring(0, 5) === '2:Del' ||
      group.mutation.substring(0, 5) === '3:Del' ||
      group.mutation.substring(0, 5) === '4:Del'
        ? group.mutation.substring(0, 1)
        : group.mutation,
    showarrow: false,
    font: { color: 'white' },
    align: 'center',
  }));

  const yTitleAnnotation = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'middle',
    yanchor: 'middle',
    align: 'center',
    x: -0.051,
    y: 0.5,
    text: '<b>Relative contribution</b>',
    textangle: -90,
    showarrow: false,
  };
  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: 700,
    title: '<b>RSS = ' + rss + '; Cosine Simularity =' + cosine + '</b>',
    grid: {
      rows: 3,
      column: 1,
    },
    autosize: true,
    xaxis: {
      showticklabels: true,
      showline: true,
      tickfont: { size: 11 },
      tickmode: 'array',
      tickvals: indelNames.map((_, i) => i),
      ticktext: indelNames.map((e) => e.index),
      tickangle: -90,
      linecolor: 'black',
      linewidth: 1,
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
      autorange: false,
      range: [0, maxMutations * 1.3],
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
      autorange: false,
      range: [0, maxMutations * 1.3],
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      domain: [0.67, 1],
    },

    shapes: [
      ...shapeTop,
      shapeRight3,
      shapeRight2,
      shapeRight1,
      ...shapeLine3,
      ...shapeLine2,
      ...shapeLine1,
    ],
    annotations: [
      ...mutationAnnotation,
      annotationLabelRight3,
      annotationLabelRight2,
      annotationLabelRight1,
      yTitleAnnotation,
      //...xLabelAnnotation,
    ],
  };

  return { traces, layout };
}
