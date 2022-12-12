import {
  groupDataByMutation,
  getMaxMutations,
} from '../mutationalProfiles/utils';

export default function pcBetweenSamples_SBS(samples, apiData) {
  const colors = {
    'C>A': '#03BCEE',
    'C>G': 'black',
    'C>T': '#E32926',
    'T>A': '#CAC9C9',
    'T>C': '#A1CE63',
    'T>G': '#EBC6C4',
  };
  const mutationRegex = /\[(.*)\]/;
  function mutationGroupSort(a, b) {
    const order = Object.keys(colors);
    return order.indexOf(a.mutation) - order.indexOf(b.mutation);
  }

  // get samples from sample array args
  const [sample1, sample2] = samples;
  // filter api data by samples
  const sampleApiData1 = apiData.filter(
    (e) => e.sample == sample1 || e.signatureName == sample1
  );
  const sampleApiData2 = apiData.filter(
    (e) => e.sample == sample2 || e.signatureName == sample2
  );

  // get total mutations per sample
  const totalMutations1 = sampleApiData1.reduce(
    (total, e) => total + (e.mutations || e.contribution),
    0
  );
  const totalMutations2 = sampleApiData2.reduce(
    (total, e) => total + (e.mutations || e.contribution),
    0
  );

  // get max mutations per sample
  const maxMutation1 = getMaxMutations(sampleApiData1) / totalMutations1;
  const maxMutation2 = getMaxMutations(sampleApiData2) / totalMutations2;
  const maxMutations = Math.max(maxMutation1, maxMutation2);

  // normalize mutations per sample
  const normalizedSample1 = sampleApiData1.map((e) => ({
    ...e,
    ...(e.mutations && {
      mutations: e.mutations / totalMutations1,
    }),
    ...(e.contribution && {
      contribution: e.contribution / totalMutations1,
    }),
  }));
  const normalizedSample2 = sampleApiData2.map((e) => ({
    ...e,
    ...(e.mutations && {
      mutations: e.mutations / totalMutations2,
    }),
    ...(e.contribution && {
      contribution: e.contribution / totalMutations2,
    }),
  }));

  const groupSamples1 = groupDataByMutation(
    normalizedSample1,
    mutationRegex,
    mutationGroupSort
  );
  const groupSamples2 = groupDataByMutation(
    normalizedSample2,
    mutationRegex,
    mutationGroupSort
  );

  // calcualte difference between samples
  const sampleDifferenceData = normalizedSample1.map((e, i) => ({
    ...e,
    ...(e.mutations && {
      mutations: e.mutations - normalizedSample2[i].mutations,
    }),
    ...(e.contribution && {
      contribution: e.contribution - normalizedSample2[i].contribution,
    }),
  }));
  const groupDifference = groupDataByMutation(
    sampleDifferenceData,
    mutationRegex,
    mutationGroupSort
  );

  const squarediff = sampleDifferenceData.map((e) =>
    Math.pow(e.mutations || e.contribution, 2)
  );
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

  const cosine = cosineSimilarity(
    normalizedSample1.map((e) => e.mutations || e.contribution),
    normalizedSample2.map((e) => e.mutations || e.contribution)
  ).toFixed(3);

  const mutationTypeNames = groupSamples1
    .map((group) =>
      group.data.map((e) => ({
        mutation: group.mutation,
        mutationType: e.mutationType,
      }))
    )
    .flat();

  const sampleTrace1 = groupSamples1.map((group, groupIndex, array) => ({
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
    y: group.data.map((e) => e.mutations || e.contribution),
    hoverinfo: 'x+y',
    showlegend: false,
    yaxis: 'y3',
  }));

  const sampleTrace2 = groupSamples2.map((group, groupIndex, array) => ({
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
    y: group.data.map((e) => e.mutations || e.contribution),
    hoverinfo: 'x+y',
    showlegend: false,
    yaxis: 'y2',
  }));
  const differenceTrace = groupDifference.map((group, groupIndex, array) => ({
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
    y: group.data.map((e) => e.mutations || e.contribution),
    hoverinfo: 'x+y',
    showlegend: false,
  }));
  const traces = [...differenceTrace, ...sampleTrace2, ...sampleTrace1];

  const shapeTop = groupSamples1.map((group, groupIndex, array) => ({
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

  const shapeLine3 = groupSamples1.map((group, groupIndex, array) => ({
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

  const shapeLine2 = groupSamples2.map((group, groupIndex, array) => ({
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

  const shapeLine1 = groupDifference.map((group, groupIndex, array) => ({
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

  const sample1Label = {
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

  const sample2Label = {
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

  const differenceLabel = {
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
  const mutationAnnotation = groupSamples1.map((group, groupIndex, array) => ({
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
    title: '<b>RSS = ' + rss + '; Cosine Simularity = ' + cosine + '</b>',
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
      sample1Label,
      sample2Label,
      differenceLabel,
      yTitleAnnotation,
    ],
  };
  return { traces, layout };
}
