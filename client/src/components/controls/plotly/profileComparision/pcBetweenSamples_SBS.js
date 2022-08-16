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
  console.log(sample1);
  console.log(sample2);

  const group1 = groupBy(sample1, 'mutationType');
  Object.keys(group1);
  console.log(group1);
  const group2 = groupBy(sample2, 'mutationType');
  Object.keys(group2);
  console.log(group2);

  let sampleDifferences = [];

  for (let mutationType of Object.keys(group1)) {
    const a = group1[mutationType][0];
    const b = group2[mutationType][0];
    const mutations = a.mutations - b.mutations;
    //const cancer = a.cancer;
    sampleDifferences.push({ mutationType, mutations });
  }
  console.log(sampleDifferences);

  const groupByMutation1 = sample1.reduce((acc, e, i) => {
    const mutationRegex = /\[(.*)\]/;
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
    const mutationRegex = /\[(.*)\]/;
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

  // const groupByMutation3 = groupBy(
  //   sampleDifferences,
  //   (s) => /\[(.*)\]/.match(s.mutationType)[1]
  // );
  const groupByMutation3 = sampleDifferences.reduce((acc, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.mutationType.match(mutationRegex)[1];
    acc[mutation] = acc[mutation] ? [...acc[mutation], e] : [e];
    return acc;
  }, {});
  console.log(groupByMutation3);

  const sample3data = Object.entries(groupByMutation3).map(
    ([mutation, data]) => ({
      mutation,
      data,
    })
  );
  console.log(sample3data);

  const totalMutations1 = sample1data.reduce(
    (total, mutation) =>
      total +
      mutation.data.reduce((mutationSum, e) => mutationSum + e.mutations, 0),
    0
  );
  const maxMutation1 = Math.max(
    ...sample1data
      .map((mutation) => mutation.data.map((e) => e.mutations))
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
      .map((mutation) => mutation.data.map((e) => e.mutations))
      .flat()
  );

  console.log(totalMutations1);
  console.log(maxMutation1);
  console.log(totalMutations2);
  console.log(maxMutation2);

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
    y: group.data.map((e) => e.mutations),
    hoverinfo: 'x+y',
    showlegend: false,
    yaxis: 'y3',
  }));
  console.log(trace1);
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
    y: group.data.map((e) => e.mutations),
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
  console.log(traces);

  const shapeTop = sample1data.map((group, groupIndex, array) => ({
    type: 'rect',
    xref: 'x',
    yref: 'paper',
    x0: array
      .slice(0, groupIndex)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.35),
    x1: array
      .slice(0, groupIndex + 1)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.65),
    y0: 1.05,
    y1: 1.01,
    fillcolor: colors[group.mutation],
    line: {
      width: 0,
    },
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
        family: 'Courier New, monospace',
        size: 8,
      },
      showgrid: true,
      gridcolor: '#F5F5F5',
      domain: [0, 0.33],
    },
    yaxis2: {
      autorange: true,
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      tickfont: {
        family: 'Courier New, monospace',
        size: 8,
      },
      domain: [0.34, 0.66],
    },
    yaxis3: {
      autorange: true,
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      tickfont: {
        family: 'Courier New, monospace',
        size: 8,
      },
      domain: [0.67, 1],
    },
    shapes: shapeTop,
  };

  return { traces, layout };
}
