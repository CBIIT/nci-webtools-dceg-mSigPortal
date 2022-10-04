export function groupDataByMutation(apiData, regex) {
  const groupByMutation = apiData.reduce((acc, e) => {
    const mutation = e.mutationType.match(regex)[1];
    acc[mutation] = acc[mutation] ? [...acc[mutation], e] : [e];
    return acc;
  }, {});

  return Object.entries(groupByMutation).map(([mutation, data]) => ({
    mutation,
    data,
  }));
}

export function getTotalMutations(apiData) {
  return apiData.reduce(
    (total, e) => total + e.mutations || e.contribution || 0,
    0
  );
}

export function getMaxMutations(apiData) {
  return Math.max(
    ...apiData
      .filter((e) => e.mutations || e.contribution)
      .map((e) => e.mutations || e.contribution)
  );
}

export function createSampleAnnotation(apiData) {
  const totalMutations = getTotalMutations(apiData);
  return {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: 0.01,
    y: 0.88,
    text: apiData[0].sample
      ? `<b>${
          apiData[0].sample
        }: ${totalMutations.toLocaleString()} Substiutions</b>`
      : `<b>${apiData[0].signatureName}</b>`,
    showarrow: false,
    font: {
      size: 24,
      family: 'Arial',
    },
    align: 'center',
  };
}

export function createMutationShapes(data, colors) {
  return data.map((group, groupIndex, array) => ({
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
}

export function createMutationAnnotations(data, appendedText = '') {
  return data.map((group, groupIndex, array) => ({
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x:
      array
        .slice(0, groupIndex)
        .reduce((lastIndex, b) => lastIndex + b.data.length, 0) +
      (group.data.length - 1) * 0.5,
    y: 1.05,
    text: `<b>${group.mutation + appendedText}</b>`,
    showarrow: false,
    font: { size: 18 },
    align: 'center',
  }));
}
