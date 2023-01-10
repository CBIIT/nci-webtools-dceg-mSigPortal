export function groupDataByMutation(
  data,
  groupRegex,
  mutationGroupSort = false,
  mutationTypeSort = false
) {
  const groupByMutation = data.reduce((acc, e) => {
    const mutation = e.mutationType.match(groupRegex)[1];
    acc[mutation] = acc[mutation] ? [...acc[mutation], e] : [e];
    return acc;
  }, {});

  const groupedData = Object.entries(groupByMutation).map(
    ([mutation, data]) => ({
      mutation,
      data: mutationTypeSort ? data.sort(mutationTypeSort) : data,
    })
  );

  return mutationGroupSort ? groupedData.sort(mutationGroupSort) : groupedData;
}

export function getTotalMutations(data) {
  return data.reduce(
    (total, e) => total + (e.mutations || e.contribution || 0),
    0
  );
}

export function getMaxMutations(data) {
  return Math.max(...data.map((e) => e.mutations || e.contribution || 0));
}

export function getRss(sampleDifferenceData) {
  const squareDiff = sampleDifferenceData.map((e) =>
    Math.pow(e.mutations || e.contribution || 0, 2)
  );
  return squareDiff.reduce((a, b, i) => a + b, 0).toExponential(3);
}

export function getCosineSimilarity(data1, data2) {
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
  return cosineSimilarity(
    data1.map((e) => e.mutations || e.contribution || 0),
    data2.map((e) => e.mutations || e.contribution || 0)
  ).toFixed(3);
}

export function compareProfiles(
  data1,
  data2,
  colors,
  mutationRegex,
  formatMutationLabels,
  formatTickLabels,
  tickAngle = -90
) {
  const sample1 = data1[0].sample || data1[0].signatureName;
  const sample2 = data2[0].sample || data2[0].signatureName;

  const mutationGroupSort = (a, b) => {
    const order = Object.keys(colors);
    return order.indexOf(a.mutation) - order.indexOf(b.mutation);
  };

  // get total mutations per sample
  const totalMutations1 = getTotalMutations(data1);
  const totalMutations2 = getTotalMutations(data2);

  // get max mutations per sample
  const maxMutation1 = getMaxMutations(data1) / totalMutations1;
  const maxMutation2 = getMaxMutations(data2) / totalMutations2;
  const maxMutations = Math.max(maxMutation1, maxMutation2);

  // normalize mutations per sample
  const normalizedSample1 = data1.map((e) => ({
    ...e,
    ...(e.mutations >= 0 && {
      mutations: e.mutations / totalMutations1,
    }),
    ...(e.contribution >= 0 && {
      contribution: e.contribution / totalMutations1,
    }),
  }));
  const normalizedSample2 = data2.map((e) => ({
    ...e,
    ...(e.mutations >= 0 && {
      mutations: e.mutations / totalMutations2,
    }),
    ...(e.contribution >= 0 && {
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
    mutations:
      (e.mutations || e.contribution || 0) -
      (normalizedSample2[i].mutations ||
        normalizedSample2[i].contribution ||
        0),
  }));
  const groupDifference = groupDataByMutation(
    sampleDifferenceData,
    mutationRegex,
    mutationGroupSort
  );

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
    y: group.data.map((e) => e.mutations || e.contribution || 0),
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
    y: group.data.map((e) => e.mutations || e.contribution || 0),
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
    y: group.data.map((e) => e.mutations || e.contribution || 0),
    hoverinfo: 'x+y',
    showlegend: false,
  }));
  const traces = [...differenceTrace, ...sampleTrace2, ...sampleTrace1];

  const rss = getRss(sampleDifferenceData);
  const cosineSimilarity = getCosineSimilarity(
    normalizedSample1,
    normalizedSample2
  );

  const tickLabels = formatTickLabels(groupSamples1);

  const mutationLabelBox = groupSamples1.map((group, groupIndex, array) => ({
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
  const sampleBorder1 = groupSamples1.map((group, groupIndex, array) => ({
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

  const sampleBorder2 = groupSamples2.map((group, groupIndex, array) => ({
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

  const differenceBorder = groupDifference.map((group, groupIndex, array) => ({
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
  const differencLabelBox1 = {
    type: 'rect',
    xref: 'paper',
    yref: 'paper',
    x0: 1,
    x1: 1.02,
    y0: 0.67,
    y1: 1,
    fillcolor: '#F0F0F0',
  };
  const sampleLabelBox2 = {
    type: 'rect',
    xref: 'paper',
    yref: 'paper',
    x0: 1,
    x1: 1.02,
    y0: 0.34,
    y1: 0.66,
    fillcolor: '#F0F0F0',
  };

  const differenceLabelBox = {
    type: 'rect',
    xref: 'paper',
    yref: 'paper',
    x0: 1,
    x1: 1.02,
    y0: 0,
    y1: 0.33,
    fillcolor: '#F0F0F0',
  };

  const sampleLabel1 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'top',
    yanchor: 'middle',
    align: 'center',
    x: 1.017,
    y: 0.835,
    text: sample1.length > 16 ? sample1.substring(0, 16) + '...' : sample1,
    textangle: 90,
    showarrow: false,
  };

  const sampleLabel2 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'top',
    yanchor: 'middle',
    align: 'center',
    x: 1.017,
    y: 0.505,
    text: sample2.length > 16 ? sample2.substring(0, 16) + '...' : sample2,
    textangle: 90,
    showarrow: false,
  };

  const differenceLabel = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'top',
    yanchor: 'middle',
    align: 'center',
    x: 1.017,
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
    y: 1.005,
    text: formatMutationLabels(group),
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
    x: -0.04,
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
    title:
      '<b>RSS = ' + rss + '; Cosine Similarity = ' + cosineSimilarity + '</b>',
    xaxis: {
      showline: true,
      tickangle: tickAngle,
      tickfont: { family: 'Courier New, monospace' },
      tickmode: 'array',
      tickvals: tickLabels.map((_, i) => i),
      ticktext: tickLabels.map((e) => e),
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
      ...mutationLabelBox,
      differencLabelBox1,
      sampleLabelBox2,
      differenceLabelBox,
      ...sampleBorder1,
      ...sampleBorder2,
      ...differenceBorder,
    ],
    annotations: [
      ...mutationAnnotation,
      sampleLabel1,
      sampleLabel2,
      differenceLabel,
      yTitleAnnotation,
    ],
  };

  return { traces, layout };
}
