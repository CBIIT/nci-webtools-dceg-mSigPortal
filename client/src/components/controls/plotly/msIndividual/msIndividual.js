import { groupBy } from 'lodash';
import { round } from '../../utils/utils';
import {
  groupDataByMutation,
  getTotalMutations,
  getMaxMutations,
  getCosineSimilarity,
  getRss,
} from '../profileComparision/profileComparison';

export function MsIndividualComparison(
  data,
  arg,
  colors,
  mutationRegex,
  formatMutationLabels,
  formatTickLabels,
  tickAngle = -90
) {
  console.log('MS Individual Plot');
  console.log(data);
  console.log(arg);
  const exposureData = data[0].data;
  console.log(exposureData);
  const signatureData = data[1].data;
  console.log(signatureData);
  const segmatrixData = data[2].data;
  console.log(segmatrixData);

  const exposure_groupBySignature = groupBy(
    exposureData.filter((o) => o['exposure'] > 0.01),
    'signatureName'
  );
  console.log(exposure_groupBySignature);

  const signatureNames = Object.keys(exposure_groupBySignature).map((e) => e);
  console.log(signatureNames);
  const exposureSum = Object.values(exposure_groupBySignature)
    .flat()
    .reduce((n, { exposure }) => n + exposure, 0);
  console.log(exposureSum);

  const percentSignature = Object.values(exposure_groupBySignature).map(
    (e) => ({
      signatureName: e[0].signatureName,
      exposure: e[0].exposure,
      exposureSum: exposureSum,
      percent: e[0].exposure / exposureSum,
    })
  );

  console.log(percentSignature);

  // let ptext = '';
  // for (let i = 0; i < percentSignature.length; i++) {
  //   ptext +=
  //     (percentSignature[i].percent * 100).toFixed(1) +
  //     '%*' +
  //     percentSignature[i].signatureName +
  //     ' + ';
  // }
  const ptext = percentSignature
    .map(
      (signature) =>
        `${(signature.percent * 100).toFixed(1)}%*${signature.signatureName} + `
    )
    .join('');

  console.log(ptext);
  const signature_groupBySignature = groupBy(
    signatureData.filter((e) => signatureNames.includes(e.signatureName)),
    'signatureName'
  );
  console.log(signature_groupBySignature);
  console.log(signatureNames.length);

  // let plotYrange;
  // if (signatureNames.length > 6) {
  //   plotYrange = 0.8;
  // } else if (signatureNames.length === 6) {
  //   plotYrange = 0.7;
  // } else if (signatureNames.length === 5) {
  //   plotYrange = 0.65;
  // } else if (signatureNames.length === 4) {
  //   plotYrange = 0.6;
  // } else {
  //   plotYrange = 0.6;
  // }
  // const divide = plotYrange / signatureNames.length;

  const plotYrange =
    signatureNames.length > 6
      ? 0.7
      : signatureNames.length === 6
      ? 0.65
      : signatureNames.length === 5
      ? 0.6
      : 0.55;

  const divide = plotYrange / signatureNames.length;
  console.log(divide);

  const signatureDataFilter = Object.values(signature_groupBySignature).flat();
  console.log(signatureDataFilter);

  const signatureDataFiltergroupBymutationTypes = groupBy(
    signatureDataFilter,
    'mutationType'
  );

  console.log(signatureDataFiltergroupBymutationTypes);

  const mutationTypes = Object.keys(
    signatureDataFiltergroupBymutationTypes
  ).map((e) => e);
  // console.log(mutationTypes);

  const seqmatrix_groupByMutationType = groupBy(
    segmatrixData.filter((e) => mutationTypes.includes(e.mutationType)),
    'mutationType'
  );
  console.log(seqmatrix_groupByMutationType);

  const seqmatrixDataFilter = Object.values(
    seqmatrix_groupByMutationType
  ).flat(); //original data for the comparison
  console.log(seqmatrixDataFilter);

  const mutationGroupSort = (a, b) => {
    const order = Object.keys(colors);
    return order.indexOf(a.mutation) - order.indexOf(b.mutation);
  };

  const totalMutationsOriginal = getTotalMutations(seqmatrixDataFilter);
  //const totalMutationsDestructed = getTotalMutations(data2);

  const maxMutationOriginal =
    getMaxMutations(seqmatrixDataFilter) / totalMutationsOriginal;

  console.log(maxMutationOriginal);
  const normalizedOriginal = seqmatrixDataFilter.map((e) => ({
    ...e,
    ...(e.mutations >= 0 && {
      mutations: e.mutations / totalMutationsOriginal,
    }),
    ...(e.contribution >= 0 && {
      contribution: e.contribution / totalMutationsOriginal,
    }),
  }));
  console.log(normalizedOriginal);

  const groupOriginal = groupDataByMutation(
    normalizedOriginal,
    mutationRegex,
    mutationGroupSort
  );
  console.log(groupOriginal);

  const arraySignatureData = Object.values(signature_groupBySignature).map(
    (e) => e
  );

  console.log(arraySignatureData);
  const arraySignatureDataFlat = arraySignatureData.flat();
  console.log(arraySignatureDataFlat);
  const destructedData = [];

  for (let i = 0; i < percentSignature.length; i++) {
    for (let j = 0; j < arraySignatureDataFlat.length; j++) {
      if (
        arraySignatureDataFlat[j].signatureName ===
        percentSignature[i].signatureName
      ) {
        let n = {
          signatureName: arraySignatureDataFlat[j].signatureName,
          mutationType: arraySignatureDataFlat[j].mutationType,
          mutations:
            arraySignatureDataFlat[j].contribution *
            percentSignature[i].percent,
        };
        destructedData.push(n);
      }
    }
  }
  // const destructedData = percentSignature
  //   .map(({ signatureName, percent }) => {
  //     const matchingData = arraySignatureDataFlat.find(
  //       ({ signatureName: name }) => name === signatureName
  //     );
  //     if (!matchingData) return null;
  //     const { mutationType, contribution } = matchingData;
  //     return {
  //       signatureName,
  //       mutationType,
  //       mutations: contribution * percent,
  //     };
  //   })
  //   .filter(Boolean);

  console.log(destructedData);

  const groupByMutationType_destructed = groupBy(
    destructedData,
    'mutationType'
  );
  console.log(Object.values(groupByMutationType_destructed).length);
  let newDestructedData = [];
  const groupByMutationType_destructed_value = Object.values(
    groupByMutationType_destructed
  );
  for (let i = 0; i < groupByMutationType_destructed_value.length; i++) {
    let n = {
      mutations: getTotalMutations(groupByMutationType_destructed_value[i]),
      mutationType: groupByMutationType_destructed_value[i][0].mutationType,
      signatureName: groupByMutationType_destructed_value[i][0].signatureName,
    };
    newDestructedData.push(n);
  }
  console.log(newDestructedData);
  // const newDestructedData = Object.values(
  //   groupBy(destructedData, 'mutationType')
  // ).map((group) => ({
  //   mutations: getTotalMutations(group),
  //   mutationType: group[0].mutationType,
  //   signatureName: group[0].signatureName,
  // }));
  // console.log(newDestructedData);

  const groupDestructed = groupDataByMutation(
    newDestructedData,
    mutationRegex,
    mutationGroupSort
  );
  console.log(groupDestructed);

  // get total mutations per sample
  const totalMutations1 = getTotalMutations(normalizedOriginal);
  const totalMutations2 = getTotalMutations(newDestructedData);
  console.log(totalMutations1);
  console.log(totalMutations2);

  // get max mutations per sample
  const maxMutation1 = getMaxMutations(normalizedOriginal) / totalMutations1;
  const maxMutation2 = getMaxMutations(newDestructedData) / totalMutations2;
  const maxMutations = Math.max(maxMutation1, maxMutation2);
  console.log(maxMutations);
  // --- Top subplots : original, destructed, different
  const sampleTraceOriginal = groupOriginal.map((group, groupIndex, array) => ({
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
    xaxis: 'x2',
    yaxis: 'y12',
  }));
  console.log(sampleTraceOriginal);
  const sampleTraceDestructed = groupDestructed.map(
    (group, groupIndex, array) => ({
      name: group.mutations,
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
      xaxis: 'x2',
      yaxis: 'y11',
    })
  );
  console.log(sampleTraceDestructed);

  const differenceTrace = sampleTraceOriginal.map((trace, traceIndex) => ({
    ...trace,
    y: trace.y.map((e, i) => e - sampleTraceDestructed[traceIndex].y[i]),
    yaxis: 'y10',
    axis: 'x2',
  }));

  const sampleBorder1 = groupOriginal.map((group, groupIndex, array) => ({
    type: 'rect',
    xref: 'x',
    yref: 'paper',
    x0: array
      .slice(0, groupIndex)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
    x1: array
      .slice(0, groupIndex + 1)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
    y0: 1 - (1 - plotYrange - 0.05) / 3,
    y1: 1,
    line: {
      width: 1,
    },
  }));

  const sampleBorder2 = groupDestructed.map((group, groupIndex, array) => ({
    type: 'rect',
    xref: 'x',
    yref: 'paper',
    x0: array
      .slice(0, groupIndex)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
    x1: array
      .slice(0, groupIndex + 1)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
    y0: 1 - ((1 - plotYrange - 0.05) / 3) * 2 - 0.01,
    y1: 1 - (1 - plotYrange - 0.05) / 3 - 0.01,

    line: {
      width: 1,
    },
  }));

  const differenceBorder = groupDestructed.map((group, groupIndex, array) => ({
    type: 'rect',
    xref: 'x',
    yref: 'paper',
    x0: array
      .slice(0, groupIndex)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
    x1: array
      .slice(0, groupIndex + 1)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
    y0: 1 - ((1 - plotYrange - 0.05) / 3) * 3 - 0.02,
    y1: 1 - ((1 - plotYrange - 0.05) / 3) * 2 - 0.02,

    line: {
      width: 1,
    },
  }));

  const sampleLabel1 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'center',
    yanchor: 'middle',
    align: 'center',
    x: 1.017,
    y: 1 - (1 - plotYrange - 0.05) / 3 + divide / 2,

    text: 'Original',
    textangle: 90,
    showarrow: false,
    width: 100,
  };

  const sampleLabel2 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'center',
    yanchor: 'middle',
    align: 'center',
    x: 1.017,
    y: 1 - ((1 - plotYrange - 0.05) / 3) * 2 - 0.01 + divide / 2,
    text: 'Reconstructed',
    textangle: 90,
    showarrow: false,
  };

  const sampleLabelDiff = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'center',
    yanchor: 'middle',
    align: 'center',
    x: 1.017,
    y: 1 - ((1 - plotYrange - 0.05) / 3) * 3 - 0.02 + divide / 2,
    text: 'Difference',
    textangle: 90,
    showarrow: false,
    height: 15,
    valign: 'top',
  };

  //-------- under subplot -----------//

  const contributionGroupSort = (a, b) => {
    const order = Object.keys(colors);
    return order.indexOf(a.contribution) - order.indexOf(b.contribution);
  };
  console.log(arraySignatureData);

  let groupSamples = [];
  for (let i = 0; i < arraySignatureData.length; i++) {
    groupSamples.push(
      groupDataByMutation(
        arraySignatureData[i],
        mutationRegex,
        contributionGroupSort
      )
    );
  }
  console.log(groupSamples);

  const tracesArray = [];
  const sampleLabels = [];
  const sampleBorders = [];
  for (let i = 0; i < groupSamples.length; i++) {
    let l = {
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      yanchor: 'middle',
      align: 'center',
      x: 1.017,
      y: divide * i + divide / 2 - 0.02,

      text: groupSamples[i][0].data[0].signatureName,
      textangle: 90,
      showarrow: false,
      width: 100,
    };

    sampleLabels.push(l);
    for (let j = 0; j < groupSamples[i].length; j++) {
      let t = {
        name: groupSamples[i][j].mutation,
        type: 'bar',
        marker: { color: colors[groupSamples[i][j].mutation] },
        x: [...groupSamples[i][j].data.keys()].map(
          (e) =>
            e +
            groupSamples[i]
              .slice(0, j)
              .reduce((lastIndex, b) => lastIndex + b.data.length, 0)
        ),
        y: groupSamples[i][j].data.map(
          (e) => e.mutations || e.contribution || 0
        ),
        hoverinfo: 'x+y',
        showlegend: false,
        yaxis: i > 0 ? 'y' + parseInt(Number(i) + Number(1)) : 'y',
      };
      tracesArray.push(t);

      let s = {
        type: 'rect',
        xref: 'x',
        yref: 'paper',
        x0: groupSamples[i]
          .slice(0, j)
          .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
        x1: groupSamples[i]
          .slice(0, j + 1)
          .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
        y0: i === 0 ? 0 : divide * i - 0.01,
        y1: divide * i + divide - 0.02,

        line: {
          width: 1,
        },
      };
      sampleBorders.push(s);
    }
  }
  console.log(sampleBorders);
  console.log(sampleTraceOriginal);
  console.log(tracesArray);
  const traces = [
    ...tracesArray,
    ...differenceTrace,
    ...sampleTraceOriginal,
    ...sampleTraceDestructed,
  ];
  console.log(traces);
  console.log(divide);

  const leftTitleAnnotation = [
    {
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      yanchor: 'middle',
      align: 'center',
      x: -0.05,
      y: plotYrange / 2,
      text: '<b>Relative contribution</b>',
      font: { size: 15 },
      textangle: -90,
      showarrow: false,
    },
    {
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      yanchor: 'middle',
      align: 'center',
      x: -0.05,
      y: plotYrange + (1 - plotYrange) / 2,
      text: '<b>Relative contribution</b>',
      font: { size: 15 },
      textangle: -90,
      showarrow: false,
    },
  ];
  const mutationAnnotation0 = groupSamples[0].map(
    (group, groupIndex, array) => ({
      xref: 'x',
      yref: 'paper',
      xanchor: 'bottom',
      yanchor: 'bottom',
      x:
        array
          .slice(0, groupIndex)
          .reduce((lastIndex, b) => lastIndex + b.data.length, 0) +
        (group.data.length - 1) * 0.5,
      y: plotYrange - 0.02,
      text: formatMutationLabels(group),
      showarrow: false,
      font: { color: 'white' },
      align: 'center',
    })
  );

  const mutationAnnotation1 = groupOriginal.map((group, groupIndex, array) => ({
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

  const percentAnnotation = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: 0.5,
    y: -0.07,
    text: '<b>Original Profile = ' + ptext.slice(0, -2) + '</b>',
    font: { size: 15 },
    showarrow: false,
    align: 'center',
  };

  const mutationLabelBox0 = groupSamples[0].map((group, groupIndex, array) => ({
    type: 'rect',
    xref: 'x',
    yref: 'paper',
    x0: array
      .slice(0, groupIndex)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
    x1: array
      .slice(0, groupIndex + 1)
      .reduce((lastIndex, e) => lastIndex + e.data.length, -0.5),
    y0: plotYrange - 0.02,
    y1: plotYrange + 0.005,
    fillcolor: colors[group.mutation],
    line: {
      width: 1,
    },
  }));

  const mutationLabelBox1 = groupOriginal.map((group, groupIndex, array) => ({
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
    y1: 1.035,
    fillcolor: colors[group.mutation],
    line: {
      width: 1,
    },
  }));
  const tickLabels = formatTickLabels(groupSamples[0]);

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: 1080,
    autosize: true,

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
      anchor: 'y',
    },
    xaxis2: {
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
      anchor: 'y10',
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
      domain: [0, divide - 0.02],
      anchor: 'x',
    },
    yaxis2: {
      autorange: true,
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      domain: [divide * 1 - 0.01, divide * 1 + divide - 0.02],
      anchor: 'x',
    },
    yaxis3: {
      autorange: true,
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      domain: [divide * 2 - 0.01, divide * 2 + divide - 0.02],
      anchor: 'x',
    },
    yaxis4: {
      autorange: true,
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      domain: [divide * 3 - 0.01, divide * 3 + divide - 0.02],
      anchor: 'x',
    },
    yaxis5: {
      autorange: true,
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      domain: [divide * 4 - 0.01, divide * 4 + divide - 0.02],
      anchor: 'x',
    },
    yaxis6: {
      autorange: true,
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      domain: [divide * 5 - 0.01, divide * 5 + divide - 0.02],
    },

    yaxis10: {
      autorange: false,
      //range: [-0.02, 0.02],
      range:
        maxMutation1 - maxMutation2 > 0
          ? [
              -1 * (maxMutation1 - maxMutation2) * 4,
              (maxMutation1 - maxMutation2) * 4,
            ]
          : [
              1 * (maxMutation1 - maxMutation2) * 4,
              -1 * (maxMutation1 - maxMutation2) * 4,
            ],
      linecolor: '#D3D3D3',
      linewidth: 1,
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      showgrid: true,
      gridcolor: '#F5F5F5',
      domain: [
        1 - ((1 - plotYrange - 0.05) / 3) * 3 - 0.015,
        1 - ((1 - plotYrange - 0.05) / 3) * 2 - 0.02,
      ],
      anchor: 'x2',
    },
    yaxis11: {
      autorange: false,
      range: [0, maxMutations * 1.2],
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      domain: [
        1 - ((1 - plotYrange - 0.05) / 3) * 2 - 0.01,
        1 - (1 - plotYrange - 0.05) / 3 - 0.01,
      ],
      anchor: 'x2',
    },
    yaxis12: {
      autorange: false,
      range: [0, maxMutations * 1.2],
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      domain: [1 - (1 - plotYrange - 0.05) / 3, 1],
      anchor: 'x2',
    },
    shapes: [
      ...mutationLabelBox0,
      ...mutationLabelBox1,
      ...sampleBorders,
      ...sampleBorder1,
      ...sampleBorder2,
      ...differenceBorder,
    ],
    annotations: [
      ...mutationAnnotation0,
      ...mutationAnnotation1,
      ...sampleLabels,
      sampleLabel1,
      sampleLabel2,
      sampleLabelDiff,
      ...leftTitleAnnotation,
      percentAnnotation,
    ],

    margin: {
      l: -0.5,
    },
  };

  console.log(layout);
  return { traces, layout };
}
