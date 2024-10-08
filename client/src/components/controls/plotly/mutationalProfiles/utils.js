export function groupDataByMutation(
  apiData,
  groupRegex,
  mutationGroupSort = false,
  mutationTypeSort = false
) {
  const groupByMutation = apiData.reduce((acc, e) => {
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

export function getTotalMutations(apiData) {
  return apiData.reduce((total, e) => {
    const mutations = Number(e.mutations) || 0;
    const contribution = Number(e.contribution) || 0;
    return total + mutations + contribution;
  }, 0);
}

export function getMaxMutations(data) {
  return Math.max(
    ...data.map((e) => Number(e.mutations) || Number(e.contribution) || 0)
  );
}

export function createSampleAnnotation(apiData, text = '', yPos = 0.88) {
  const totalMutations = getTotalMutations(apiData);
  return {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: 0.01,
    y: yPos,
    text:
      apiData[0].sample && parseFloat(totalMutations).toFixed(2) > 1
        ? `<b>${apiData[0].sample}: ${totalMutations.toLocaleString()} ${
            text || apiData[0].profile == 'ID' ? 'Indels' : 'Substitutions'
          }</b>`
        : apiData[0].sample && totalMutations <= 1.1
        ? `<b>${apiData[0].sample}</b>`
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
  function getAnnotationFontSize() {
    if (window.innerWidth < 768) {
      return 12;
    } else if (window.innerWidth >= 768 && window.innerWidth < 1024) {
      return 14;
    } else {
      return 18; // Default font size for larger screens
    }
  }

  const fontSize = getAnnotationFontSize();

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
    //font: { size: 18 },
    font: { size: fontSize }, // Use the dynamically determined font size
    align: 'center',
  }));
}
