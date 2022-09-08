import { groupBy } from 'lodash';

export default function pcBetweenSamples(samples, sample1, sample2, tab) {
  const colors = {
    AC: '#09BCED',
    AT: '#0266CA',
    CC: '#9FCE62',
    CG: '#006501',
    CT: '#FF9898',
    GC: '#E22925',
    TA: '#FEB065',
    TC: '#FD8000',
    TG: '#CB98FD',
    TT: '#4C0299',
  };

  const dbsdata = ['AC', 'AT', 'CC', 'CG', 'CT', 'GC', 'TA', 'TC', 'TG', 'TT'];

  const mutationRegex = /^(.{2})/;
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
  sample1data.filter((e) => dbsdata.includes(e.mutation));
  sample2data.filter((e) => dbsdata.includes(e.mutation));

  const mutationTypeNames = sample1data
    .map((group) =>
      group.data.map((e) => ({
        mutation: group.mutation,
        mutationType: e.mutationType,
      }))
    )
    .flat();

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
      mutation.data.reduce(
        (mutationSum, e) =>
          tab === 'samples'
            ? mutationSum + e.mutations
            : mutationSum + e.contribution,
        0
      ),
    0
  );
  const maxMutation2 = Math.max(
    ...sample2data
      .map((mutation) =>
        mutation.data.map((e) =>
          tab === 'samples'
            ? e.mutations / totalMutations2
            : e.contribution / totalMutations2
        )
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
    let mutations;
    tab === 'samples'
      ? (mutations =
          a.mutations / totalMutations1 - b.mutations / totalMutations2)
      : (mutations =
          a.mutations / totalMutations1 - b.contribution / totalMutations2);
    sampleDifferences.push({ mutationType, mutations });
    s1mutations.push(a.mutations / totalMutations1);
    s2mutations.push(
      tab === 'samples'
        ? b.mutations / totalMutations2
        : b.contribution / totalMutations2
    );
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

  const trace1 = sample1data.map((group, groupIndex, array) => ({
    name: group.mutation,
    type: 'bar',
    marker: { color: colors[group.mutation] },
    x: [...group.data.keys()].map(
      (e) =>
        e +
        array
          .slice(0, groupIndex)
          .reduce((lastIndex, b) => lastIndex + b.data.length, 0)
    ),
    y: group.data.map((e) => e.mutations / totalMutations1),
    hoverinfo: 'x+y',
    showlegend: false,
    yaxis: 'y3',
  }));

  const trace2 = sample2data.map((group, groupIndex, array) => ({
    name: group.mutation,
    type: 'bar',
    marker: { color: colors[group.mutation] },
    x: [...group.data.keys()].map(
      (e) =>
        e +
        array
          .slice(0, groupIndex)
          .reduce((lastIndex, b) => lastIndex + b.data.length, 0)
    ),
    y: group.data.map((e) =>
      tab === 'samples'
        ? e.mutations / totalMutations2
        : e.contribution / totalMutations2
    ),
    hoverinfo: 'x+y',
    showlegend: false,
    yaxis: 'y2',
  }));

  const trace3 = sample3data.map((group, groupIndex, array) => ({
    name: group.mutation,
    type: 'bar',
    marker: { color: colors[group.mutation] },
    x: [...group.data.keys()].map(
      (e) =>
        e +
        array
          .slice(0, groupIndex)
          .reduce((lastIndex, b) => lastIndex + b.data.length, 0)
    ),
    y: group.data.map((e) => e.mutations),
    hoverinfo: 'x+y',
    showlegend: false,
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
    fillcolor: colors[group.mutation],
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
    text:
      samples[0].length > 16 ? samples[0].substring(0, 16) + '...' : samples[0],
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
    text:
      samples[1].length > 16 ? samples[1].substring(0, 16) + '...' : samples[1],
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
    text: `<b>${group.mutation}>NN</b>`,
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
    height: 700,
    hoverlabel: { bgcolor: '#FFF' },
    grid: {
      rows: 3,
      column: 1,
    },
    title: '<b>RSS = ' + rss + '; Cosine Simularity =' + cosine + '</b>',
    xaxis: {
      showline: true,
      tickangle: -90,
      tickfont: {
        family: 'Courier New, monospace',
        color: '#A0A0A0',
      },
      tickmode: 'array',
      tickvals: mutationTypeNames.map((_, i) => i),
      ticktext: mutationTypeNames.map((e) => e.mutationType.slice(-2)),
      linecolor: '#E0E0E0',
      linewidth: 1,
      mirror: 'all',
      ticks: '',
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
      range: [0, maxMutations * 1.2],
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
      range: [0, maxMutations * 1.2],

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
    ],
  };

  return { traces, layout };
}
