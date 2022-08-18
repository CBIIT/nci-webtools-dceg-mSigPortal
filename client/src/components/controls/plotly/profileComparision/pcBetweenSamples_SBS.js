import { groupBy } from 'lodash';

export default function pcBetweenSamples(rawData, args) {
  const colors = {
    'C>A': '#03BCEE',
    'C>G': 'black',
    'C>T': '#E32926',
    'T>A': '#CAC9C9',
    'T>C': '#A1CE63',
    'T>G': '#EBC6C4',
  };
  const samples = args.sample.split(',');
  const groupBySample = groupBy(rawData, 'sample');
  const sample1 = groupBySample[samples[0]].flat();
  const sample2 = groupBySample[samples[1]].flat();

  const totalMutations1_5 = sample1.reduce(
    (a, b, i) => (i <= 5 ? a + b.mutations : a),
    0
  );
  console.log(totalMutations1_5);

  const sample1_5 = sample1.map((e, i) =>
    i <= 5 ? e.mutations / totalMutations1_5 : 0
  );
  console.log(sample1_5);

  const totalMutations2_5 = sample2.reduce(
    (a, b, i) => (i <= 5 ? a + b.mutations : a),
    0
  );
  console.log(totalMutations2_5);

  const sample2_5 = sample2.map((e, i) =>
    i <= 5 ? e.mutations / totalMutations2_5 : 0
  );
  console.log(sample2_5);

  const diffNumber = (arr1, arr2) =>
    arr1.map(function (num, idx) {
      return num - arr2[idx];
    });

  const diff = diffNumber(sample1_5, sample2_5);
  const sum_diff = diff.reduce((a, b, i) => (i <= 5 ? a + b : a), 0);
  console.log(sum_diff);
  const rss = sum_diff * sum_diff;

  const groupByMutation1 = sample1.reduce((acc, e, i) => {
    const mutationRegex = /\[(.*)\]/;
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
    const mutationRegex = /\[(.*)\]/;
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

  for (let mutationType of Object.keys(group1)) {
    const a = group1[mutationType][0];
    const b = group2[mutationType][0];
    const mutations =
      a.mutations / totalMutations1 - b.mutations / totalMutations2;
    //const cancer = a.cancer;
    sampleDifferences.push({ mutationType, mutations });
  }

  const groupByMutation3 = groupBy(
    sampleDifferences,
    (s) => s.mutationType.match(/\[(.*)\]/)[1]
  );

  const sample3data = Object.entries(groupByMutation3).map(
    ([mutation, data]) => ({
      mutation,
      data,
    })
  );

  const mutationTypeNames = sample1data
    .map((group) =>
      group.data.map((e) => ({
        mutation: group.mutation,
        mutationType: e.mutationType,
      }))
    )
    .flat();
  console.log(mutationTypeNames);
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
    y: group.data.map((e) => e.mutations / totalMutations2),
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
  console.log(trace3);
  const traces = [...trace2, ...trace3, ...trace1];

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
  };
  const shapeRight2 = {
    type: 'rect',
    xref: 'paper',
    yref: 'paper',
    x0: 1,
    x1: 1.02,
    y0: 0.34,
    y1: 0.66,
  };

  const shapeRight1 = {
    type: 'rect',
    xref: 'paper',
    yref: 'paper',
    x0: 1,
    x1: 1.02,
    y0: 0,
    y1: 0.33,
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
  console.log(annotationLabelRight3);

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
    text: `<b>${group.mutation}</b>`,
    showarrow: false,
    font: { size: 16, color: 'white' },
    align: 'center',
  }));
  function formatTickLabel(mutation, mutationType) {
    const color = colors[mutation];
    const regex = /^(.)\[(.).{2}\](.)$/;
    const match = mutationType.match(regex);
    return `${match[1]}<span style="color:${color}"><b>${match[2]}</b></span>${match[3]}`;
  }
  const layout = {
    height: 700,
    hoverlabel: { bgcolor: '#FFF' },
    grid: {
      rows: 3,
      column: 1,
    },
    xaxis: {
      showline: true,
      tickangle: -90,
      tickfont: {
        family: 'Courier New, monospace',
        color: '#A0A0A0',
      },
      tickmode: 'array',
      tickvals: mutationTypeNames.map((_, i) => i),
      ticktext: mutationTypeNames.map((e) =>
        formatTickLabel(e.mutation, e.mutationType)
      ),
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
    ],
  };

  return { traces, layout };
}
