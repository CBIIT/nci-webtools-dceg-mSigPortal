import { groupBy } from 'lodash';

export default function pcBetweenSamples(rawData, args) {
  console.log(rawData);
  console.log(args);

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

  const samples = args.sample.split(',');
  const groupBySample = groupBy(rawData, 'sample');
  const sample1 = groupBySample[samples[0]].flat();
  const sample2 = groupBySample[samples[1]].flat();
  console.log(sample1);
  console.log(sample2);
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
  console.log(sample1data);
  console.log(sample2data);
  const mutationTypeNames = sample1data
    .map((group) =>
      group.data.map((e) => ({
        mutation: group.mutation,
        mutationType: e.mutationType,
      }))
    )
    .flat();

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

  const traces = [...trace1, ...trace2];

  const layout = {
    height: 700,
    hoverlabel: { bgcolor: '#FFF' },
    grid: {
      rows: 3,
      column: 1,
    },
    //title: '<b>RSS = ' + rss + '; Cosine Simularity =' + cosine + '</b>',
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
    // shapes: [
    //   ...shapeTop,
    //   shapeRight3,
    //   shapeRight2,
    //   shapeRight1,
    //   ...shapeLine3,
    //   ...shapeLine2,
    //   ...shapeLine1,
    // ],
    // annotations: [
    //   ...mutationAnnotation,
    //   annotationLabelRight3,
    //   annotationLabelRight2,
    //   annotationLabelRight1,
    //   yTitleAnnotation,
    // ],
  };

  return { traces, layout };
}
