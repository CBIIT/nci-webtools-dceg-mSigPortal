import { groupBy } from 'lodash';
import { round, arrayContainsTerms } from '../../utils/utils';
import { colorPallet, colorPallet1 } from '../../utils/colors';

import {
  groupDataByMutation,
  getTotalMutations,
  getMaxMutations,
  getCosineSimilarity,
  getRss,
  findMaxAbsoluteYValue,
} from '../profileComparison/profileComparison';

export function MsIndividualComparison(
  data,
  arg,
  colors,
  mutationRegex,
  formatMutationLabels,
  formatTickLabels,
  tickAngle = -90
) {
  const exposureData = data[0].data;
  const signatureData = data[1].data;
  const segmatrixData = data[2].data;

  const exposure_groupBySignature = groupBy(
    exposureData.filter((o) => o['exposure'] > 0.01),
    'signatureName'
  );

  const signatureNames = Object.keys(exposure_groupBySignature).map((e) => e);
  // find the longest label to calculate extra height margin
  const longest = signatureNames.reduce(
    (a, e) => (a > e.length ? a : e.length),
    0
  );
  const extraMargin = longest < 7 ? 200 : longest * 12.5;

  const searchTerms = ['SBS'];
  const containsTerm = arrayContainsTerms(signatureNames, searchTerms);
  let signatureColors;
  containsTerm
    ? (signatureColors = colorPallet)
    : (signatureColors = colorPallet1);
  const exposureSum = Object.values(exposure_groupBySignature)
    .flat()
    .reduce((n, { exposure }) => n + exposure, 0);

  const percentSignature = Object.values(exposure_groupBySignature).map(
    (e) => ({
      signatureName: e[0].signatureName,
      exposure: e[0].exposure,
      exposureSum: exposureSum,
      percent: e[0].exposure / exposureSum,
    })
  );

  const ptext = percentSignature
    .map(
      (signature) =>
        `${(signature.percent * 100).toFixed(1)}%*${signature.signatureName} + `
    )
    .join('');

  const signature_groupBySignature = groupBy(
    signatureData.filter((e) => signatureNames.includes(e.signatureName)),
    'signatureName'
  );

  const plotYrange2 =
    signatureNames.length > 6
      ? 0.68
      : signatureNames.length === 6
      ? 0.65
      : signatureNames.length === 5
      ? 0.6
      : signatureNames.length === 4
      ? 0.55
      : signatureNames.length === 3
      ? 0.5
      : signatureNames.length === 2
      ? 0.4
      : signatureNames.length === 1
      ? 0.2
      : 0.1;
  const plotYrange1 = 1 - plotYrange2 - 0.06;
  const divide2 = plotYrange2 / signatureNames.length;
  const divide1 = plotYrange1 / 3;

  const signatureDataFiltergroupBymutationTypes = groupBy(
    Object.values(signature_groupBySignature).flat(),
    'mutationType'
  );

  const seqmatrix_groupByMutationType = groupBy(
    segmatrixData.filter((e) =>
      Object.keys(signatureDataFiltergroupBymutationTypes)
        .map((m) => m)
        .includes(e.mutationType)
    ),
    'mutationType'
  );

  const seqmatrixDataFilter = Object.values(
    seqmatrix_groupByMutationType
  ).flat(); //original data for the comparison

  const mutationGroupSort = (a, b) => {
    const order = Object.keys(colors);
    return order.indexOf(a.mutation) - order.indexOf(b.mutation);
  };

  const totalMutationsOriginal = getTotalMutations(seqmatrixDataFilter);

  const normalizedOriginal = seqmatrixDataFilter.map((e) => ({
    ...e,
    ...(e.mutations >= 0 && {
      mutations: e.mutations / totalMutationsOriginal,
    }),
    ...(e.contribution >= 0 && {
      contribution: e.contribution / totalMutationsOriginal,
    }),
  }));

  const groupOriginal = groupDataByMutation(
    normalizedOriginal,
    mutationRegex,
    mutationGroupSort
  );

  const arraySignatureData = Object.values(signature_groupBySignature).map(
    (e) => e
  );

  // const arraySignatureDataFlat = arraySignatureData.flat();
  // const destructedData = [];
  // for (let i = 0; i < percentSignature.length; i++) {
  //   for (let j = 0; j < arraySignatureDataFlat.length; j++) {
  //     if (
  //       arraySignatureDataFlat[j].signatureName ===
  //       percentSignature[i].signatureName
  //     ) {
  //       let n = {
  //         signatureName: arraySignatureDataFlat[j].signatureName,
  //         mutationType: arraySignatureDataFlat[j].mutationType,
  //         mutations:
  //           arraySignatureDataFlat[j].contribution *
  //           percentSignature[i].percent,
  //       };
  //       destructedData.push(n);
  //     }
  //   }
  // }

  const destructedData = arraySignatureData
    .flat()
    .filter((data) =>
      percentSignature.some((p) => p.signatureName === data.signatureName)
    )
    .map((data) => ({
      signatureName: data.signatureName,
      mutationType: data.mutationType,
      mutations:
        data.contribution *
        percentSignature.find((p) => p.signatureName === data.signatureName)
          .percent,
    }));

  // const groupByMutationType_destructed = groupBy(
  //   destructedData,
  //   'mutationType'
  // );
  // let newDestructedData = [];
  // const groupByMutationType_destructed_value = Object.values(
  //   groupByMutationType_destructed
  // );
  // for (let i = 0; i < groupByMutationType_destructed_value.length; i++) {
  //   let n = {
  //     mutations: getTotalMutations(groupByMutationType_destructed_value[i]),
  //     mutationType: groupByMutationType_destructed_value[i][0].mutationType,
  //     signatureName: groupByMutationType_destructed_value[i][0].signatureName,
  //   };
  //   newDestructedData.push(n);
  // }
  const newDestructedData = Object.values(
    destructedData.reduce((acc, curr) => {
      const key = curr.mutationType;
      if (!acc[key]) {
        acc[key] = {
          mutations: 0,
          mutationType: curr.mutationType,
          signatureName: curr.signatureName,
        };
      }
      acc[key].mutations += curr.mutations;
      return acc;
    }, {})
  );

  const groupDestructed = groupDataByMutation(
    newDestructedData,
    mutationRegex,
    mutationGroupSort
  );

  // get total mutations per sample
  const totalMutations1 = getTotalMutations(normalizedOriginal);
  const totalMutations2 = getTotalMutations(newDestructedData);

  // get max mutations per sample
  const maxMutation1 = getMaxMutations(normalizedOriginal) / totalMutations1;
  const maxMutation2 = getMaxMutations(newDestructedData) / totalMutations2;
  const maxMutations = Math.max(maxMutation1, maxMutation2);
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
  const differenceTraceMaxYValue = findMaxAbsoluteYValue(differenceTrace);
  const sample1Data = sampleTraceOriginal.reduce(
    (array, trace) => [...array, ...trace.y],
    []
  );
  const sample2Data = sampleTraceDestructed.reduce(
    (array, trace) => [...array, ...trace.y],
    []
  );

  const sampleDifferenceData = differenceTrace.reduce(
    (array, trace) => [...array, ...trace.y],
    []
  );
  const rss = getRss(sampleDifferenceData);
  const cosineSimilarity = getCosineSimilarity(sample1Data, sample2Data);

  //-------- under subplot -----------//

  const contributionGroupSort = (a, b) => {
    const order = Object.keys(colors);
    return order.indexOf(a.contribution) - order.indexOf(b.contribution);
  };

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
  groupSamples.reverse(); //make the lower subplot has same order as in stage

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
      y: divide2 * i + divide2 / 2 - 0.02,

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
        y0: i === 0 ? 0 : divide2 * i - 0.01,
        y1: divide2 * i + divide2 - 0.02,

        line: {
          width: 1,
        },
      };
      sampleBorders.push(s);
    }
  }

  const traces = [
    ...tracesArray,
    ...differenceTrace,
    ...sampleTraceOriginal,
    ...sampleTraceDestructed,
  ];
  // ----- Shapes -------//
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
    y0: 1 - divide1 - 0.01,
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
    y0: 1 - divide1 * 2 - 0.01,
    y1: 1 - divide1 - 0.02,

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
    y0: 1 - divide1 * 3 - 0.01,
    y1: 1 - divide1 * 2 - 0.02,

    line: {
      width: 1,
    },
  }));

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
    y0: plotYrange2 - 0.02,
    y1: plotYrange2 + 0.005,
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
    y1: 1.025,
    fillcolor: colors[group.mutation],
    line: {
      width: 1,
    },
  }));

  const sortArr = percentSignature
    .slice()
    .sort((a, b) => b.percent - a.percent);
  const percents = sortArr.map((obj) => obj.percent); // extract percent values
  const scaledPercents = percents.map((p, i) =>
    percents.slice(0, i + 1).reduce((acc, val) => acc + val)
  ); // compute cumulative sum

  const signaturePercentBox = scaledPercents.map((val, i, arr) => ({
    type: 'rect',
    xref: 'paper',
    yref: 'paper',
    y0: i === 0 ? 0 : arr[i - 1],
    y1: val,
    x0: -0.105,
    x1: -0.08,
    signatureName: sortArr[i] ? sortArr[i].signatureName : '',
    fillcolor: containsTerm
      ? signatureColors[
          sortArr[i].signatureName.replace(/^\D*/, '').replace(')', '')
        ]
      : signatureColors[i],
    line: {
      width: 0,
    },
  }));

  const signaturePercentLine = scaledPercents.map((val, i, arr) => ({
    type: 'line',
    xref: 'paper',
    yref: 'paper',
    y0: i === 0 ? val / 2 : (val - arr[i - 1]) / 2 + arr[i - 1],
    y1: i === 0 ? val / 2 : (val - arr[i - 1]) / 2 + arr[i - 1],
    x0: -0.105,
    x1: -0.125,
    signatureName: sortArr[i] ? sortArr[i].signatureName : '',

    line: {
      width: 1,
      color: containsTerm
        ? signatureColors[
            sortArr[i].signatureName.replace(/^\D*/, '').replace(')', '')
          ]
        : signatureColors[i],
    },
  }));

  //------ Annotations -------//
  const topSubplotAnnotations = [
    {
      //Sample Original
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      yanchor: 'middle',
      align: 'center',
      x: 1.017,
      y: 1 - divide1 / 2,

      text: 'Original',
      textangle: 90,
      showarrow: false,
      width: 100,
    },
    {
      //Sample Destructed
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      yanchor: 'middle',
      align: 'center',
      x: 1.017,
      y: 1 - divide1 * 1.5 - 0.015,
      text: 'Deconstructed',
      textangle: 90,
      showarrow: false,
    },
    {
      //Sample difference
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      yanchor: 'middle',
      align: 'center',
      x: 1.017,
      y: 1 - divide1 * 2.5 - 0.015,
      text: 'Difference',
      textangle: 90,
      showarrow: false,
      valign: 'top',
    },
  ];

  const titleAnnotations = [
    {
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      yanchor: 'middle',
      align: 'center',
      x: -0.05,
      y: 1 - divide1 * 1.5 - 0.015,
      text: '<b>Relative contribution</b>',
      font: { size: 16, family: 'Times New Roman' },
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
      y: plotYrange2 / 2,
      text: '<b>Relative contribution</b>',
      font: { size: 16, family: 'Times New Roman' },
      textangle: -90,
      showarrow: false,
    },
    {
      xref: 'paper',
      yref: 'paper',
      xanchor: 'bottom',
      yanchor: 'bottom',
      x: 0.5,
      y: -0.07,
      text: '<b>Original Profile = ' + ptext.slice(0, -2) + '</b>',
      font: { size: 13 },
      showarrow: false,
      align: 'center',
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
      y: plotYrange2 - 0.0185,
      text: formatMutationLabels(group),
      showarrow: false,
      font: { color: 'white', size: 13, family: 'Times New Roman' },
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
    y: 1.002,
    text: formatMutationLabels(group),
    showarrow: false,
    font: { color: 'white', size: 13, family: 'Times New Roman' },
    align: 'center',
  }));

  const signaturePercentAnnotation = scaledPercents.map((val, i, arr) => ({
    xref: 'paper',
    yref: 'paper',
    xanchor: 'center',
    yanchor: 'middle',
    val: val,
    val1: arr[i - 1],
    y: i === 0 ? val / 2 : (val - arr[i - 1]) / 2 + arr[i - 1],
    x: longest < 7 ? -0.15 : -0.2,
    signatureName: sortArr[i] ? sortArr[i].signatureName : '',
    font: {
      color: containsTerm
        ? signatureColors[
            sortArr[i].signatureName.replace(/^\D*/, '').replace(')', '')
          ]
        : signatureColors[i],
    },
    text: sortArr[i].signatureName,
    showarrow: false,

    align: 'center',
  }));
  const tickLabels = formatTickLabels(groupSamples[0]);

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: 1080,
    autosize: true,
    title: {
      text:
        '<b>Mutational Signature Association</b><br><b>RSS = ' +
        rss +
        '; Cosine Similarity = ' +
        cosineSimilarity +
        '</b>',
      font: {
        family: 'Times New Roman',
        size: 20,
      },
    },

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
      domain: [0, divide2 - 0.02],
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
      domain: [divide2 * 1 - 0.01, divide2 * 1 + divide2 - 0.02],
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
      domain: [divide2 * 2 - 0.01, divide2 * 2 + divide2 - 0.02],
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
      domain: [divide2 * 3 - 0.01, divide2 * 3 + divide2 - 0.02],
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
      domain: [divide2 * 4 - 0.01, divide2 * 4 + divide2 - 0.02],
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
      domain: [divide2 * 5 - 0.01, divide2 * 5 + divide2 - 0.02],
    },
    yaxis7: {
      autorange: true,
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      domain: [divide2 * 6 - 0.01, divide2 * 6 + divide2 - 0.02],
    },
    yaxis8: {
      autorange: true,
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      domain: [divide2 * 8 - 0.01, divide2 * 8 + divide2 - 0.02],
    },
    yaxis10: {
      autorange: false,

      range: [
        -1 * differenceTraceMaxYValue * 1.5,
        differenceTraceMaxYValue * 1.5,
      ],

      linecolor: '#D3D3D3',
      linewidth: 1,
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      showgrid: true,
      gridcolor: '#F5F5F5',

      domain: [1 - divide1 * 3 - 0.01, 1 - divide1 * 2 - 0.02],
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
      //title: { text: '<b>Relative contribution</b>' },
      domain: [1 - divide1 * 2 - 0.01, 1 - divide1 - 0.02],
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
      domain: [1 - divide1 - 0.01, 1],
      anchor: 'x2',
    },
    shapes: [
      ...mutationLabelBox0,
      ...mutationLabelBox1,
      ...sampleBorders,
      ...sampleBorder1,
      ...sampleBorder2,
      ...differenceBorder,
      ...signaturePercentBox,
      ...signaturePercentLine,
    ],
    annotations: [
      ...mutationAnnotation0,
      ...mutationAnnotation1,
      ...sampleLabels,
      ...topSubplotAnnotations,
      ...titleAnnotations,
      ...signaturePercentAnnotation,
    ],

    margin: {
      l: extraMargin,
      t: 150,
    },
  };

  return { traces, layout };
}
