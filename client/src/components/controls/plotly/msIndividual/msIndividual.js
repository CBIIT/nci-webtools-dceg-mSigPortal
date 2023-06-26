import { groupBy } from 'lodash';
import { round, arrayContainsTerms } from '../../utils/utils.js';
import { colorPallet, colorPallet1 } from '../../utils/colors.js';

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
    exposureData?.filter((o) => o['exposure'] > 0.01),
    'signatureName'
  );

  const signatureNames = Object.keys(exposure_groupBySignature).map((e) => e);

  // find the longest label to calculate extra height margin
  const longest = signatureNames.reduce(
    (a, e) => (a > e.length ? a : e.length),
    0
  );
  const extraMargin = longest < 7 ? 200 : longest * 12.5;

  //console.log(signatureNames);

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
    signatureData?.filter((e) => signatureNames.includes(e.signatureName)),
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
    segmatrixData?.filter((e) =>
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
  // console.log('maxMutation1', maxMutation1);
  // console.log('maxMutation2 ', maxMutation2);
  // console.log(maxMutations);
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
    yaxis: 'y22',
  }));

  const sampleTraceDestructed = groupDestructed.map(
    (group, groupIndex, array) => ({
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
      yaxis: 'y21',
    })
  );
  const differenceTrace = sampleTraceOriginal.map((trace, traceIndex) => ({
    ...trace,
    y: trace.y.map((e, i) => e - sampleTraceDestructed[traceIndex].y[i]),
    yaxis: 'y20',
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

  const groupSamples = [];
  for (let i = 0; i < arraySignatureData.length; i++) {
    groupSamples.push(
      groupDataByMutation(
        arraySignatureData[i],
        mutationRegex,
        mutationGroupSort
      )
    );
  }

  // const flattenedSubplotsGroup = groupSamples.map((subArray) =>
  //   subArray.flatMap((obj) => Object.values(obj.data))
  // );

  // console.log('flattenedTestGroup:', flattenedSubplotsGroup);

  // const contributionGroupSort = (a, b) => {
  //   const order = Object.keys(colors);
  //   return order.indexOf(a.contribution) - order.indexOf(b.contribution);
  // };

  // let groupSamples = [];
  // for (let i = 0; i < arraySignatureData.length; i++) {
  //   groupSamples.push(
  //     groupDataByMutation(
  //       arraySignatureData[i],
  //       mutationRegex,
  //       contributionGroupSort
  //     )
  //   );
  // }
  // groupSamples.reverse(); //make the lower subplot has same order as in stage

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
  const sampleBorder1 = groupDestructed.map((group, groupIndex, array) => ({
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

  const mutationLabelBox0 = groupSamples[0]?.map(
    (group, groupIndex, array) => ({
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
    })
  );

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

  // left side percent rectangle shape
  const signaturePercentBox = scaledPercents.map((val, i, arr) => ({
    type: 'rect',
    xref: 'paper',
    yref: 'paper',
    xanchor: 'right',
    y0: i === 0 ? 0 : arr[i - 1],
    y1: val,
    x0: -0.045,
    x1: -0.055,
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

  // left side line for percent rectangle shape
  const signaturePercentLine = scaledPercents.map((val, i, arr) => ({
    type: 'line',
    xref: 'paper',
    yref: 'paper',
    xanchor: 'right',
    y0: i === 0 ? val / 2 : (val - arr[i - 1]) / 2 + arr[i - 1],
    y1: i === 0 ? val / 2 : (val - arr[i - 1]) / 2 + arr[i - 1],
    x0: -0.055,
    x1: -0.058,
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
      xanchor: 'right',
      yanchor: 'middle',
      align: 'center',
      x: -0.025,
      y: 1 - divide1 * 1.5 - 0.015,
      text: '<b>Relative contribution</b>',
      font: { size: 16, family: 'Times New Roman' },
      textangle: -90,
      showarrow: false,
    },
    {
      xref: 'paper',
      yref: 'paper',
      xanchor: 'right',
      yanchor: 'middle',
      align: 'center',
      x: -0.025,
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
  const mutationAnnotation0 = groupSamples[0]?.map(
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

  // left percent annotation
  const signaturePercentAnnotation = scaledPercents.map((val, i, arr) => ({
    xref: 'paper',
    yref: 'paper',
    xanchor: 'right',
    yanchor: 'middle',
    val: val,
    val1: arr[i - 1],
    y: i === 0 ? val / 2 : (val - arr[i - 1]) / 2 + arr[i - 1],
    //x: longest < 7 ? -0.09 : -0.1,
    x: -0.058,
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
  const tickLabels =
    groupSamples.length > 0 ? formatTickLabels(groupSamples[0]) : '';

  const layout =
    traces.length > 0
      ? {
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
            tickvals: tickLabels.length ? tickLabels.map((_, i) => i) : '',
            ticktext: tickLabels.length ? tickLabels.map((e) => e) : '',
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
            tickvals: tickLabels.length ? tickLabels.map((_, i) => i) : '',
            ticktext: tickLabels.length ? tickLabels.map((e) => e) : '',
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
            domain: [divide2 * 7 - 0.01, divide2 * 7 + divide2 - 0.02],
          },
          yaxis9: {
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
            autorange: true,
            linecolor: '#D3D3D3',
            linewidth: 1,
            ticks: '',
            mirror: 'all',
            tickfont: {
              family: 'Arial',
            },
            domain: [divide2 * 9 - 0.01, divide2 * 9 + divide2 - 0.02],
          },
          yaxis11: {
            autorange: true,
            linecolor: '#D3D3D3',
            linewidth: 1,
            ticks: '',
            mirror: 'all',
            tickfont: {
              family: 'Arial',
            },
            domain: [divide2 * 10 - 0.01, divide2 * 10 + divide2 - 0.02],
          },
          yaxis12: {
            autorange: true,
            linecolor: '#D3D3D3',
            linewidth: 1,
            ticks: '',
            mirror: 'all',
            tickfont: {
              family: 'Arial',
            },
            domain: [divide2 * 11 - 0.01, divide2 * 11 + divide2 - 0.02],
          },
          yaxis13: {
            autorange: true,
            linecolor: '#D3D3D3',
            linewidth: 1,
            ticks: '',
            mirror: 'all',
            tickfont: {
              family: 'Arial',
            },
            domain: [divide2 * 12 - 0.01, divide2 * 12 + divide2 - 0.02],
          },
          yaxis14: {
            autorange: true,
            linecolor: '#D3D3D3',
            linewidth: 1,
            ticks: '',
            mirror: 'all',
            tickfont: {
              family: 'Arial',
            },
            domain: [divide2 * 13 - 0.01, divide2 * 13 + divide2 - 0.02],
          },

          yaxis20: {
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
          yaxis21: {
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
          yaxis22: {
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
          shapes: mutationLabelBox0
            ? [
                ...mutationLabelBox0,
                ...mutationLabelBox1,
                ...sampleBorders,
                ...sampleBorder1,
                ...sampleBorder2,
                ...differenceBorder,
                ...signaturePercentBox,
                ...signaturePercentLine,
              ]
            : [
                [
                  ...mutationLabelBox1,
                  ...sampleBorders,
                  ...sampleBorder1,
                  ...sampleBorder2,
                  ...differenceBorder,
                  ...signaturePercentBox,
                  ...signaturePercentLine,
                ],
              ],
          annotations: mutationAnnotation0
            ? [
                ...mutationAnnotation0,
                ...mutationAnnotation1,
                ...sampleLabels,
                ...topSubplotAnnotations,
                ...titleAnnotations,
                ...signaturePercentAnnotation,
              ]
            : [
                [
                  ...mutationAnnotation1,
                  ...sampleLabels,
                  ...topSubplotAnnotations,
                  ...titleAnnotations,
                  ...signaturePercentAnnotation,
                ],
              ],

          margin: {
            l: extraMargin,
            t: 150,
          },
        }
      : {
          hoverlabel: { bgcolor: '#FFF' },
          height: 250,
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
            autorange: true,
            showgrid: false,
            zeroline: false,
            showline: false,
            autotick: true,
            ticks: '',
            showticklabels: false,
          },
          yaxis: {
            autorange: true,
            showgrid: false,
            zeroline: false,
            showline: false,
            autotick: true,
            ticks: '',
            showticklabels: false,
          },
          annotations: [
            {
              text: 'No data available',
              xref: 'paper',
              yref: 'paper',
              x: 0.5,
              y: 0.5,
              showarrow: false,
              font: {
                size: 28,
                color: 'grey',
              },
            },
          ],
        };

  return { traces, layout };
}
