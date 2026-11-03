import { groupBy } from 'lodash';
import { arrayContainsTerms } from '@/components/controls/utils/utils';
import { colorPallet, colorPallet1 } from '@/components/controls/utils/colors';

import {
  groupDataByMutation,
  groupFilterDataByMutation,
  getTotalMutations,
  getMaxMutations,
  getCosineSimilarity,
  getRss,
  findMaxAbsoluteYValue,
} from '@/components/controls/plotly/profileComparison/profileComparison';

// pixel sizes; the figure height grows with the number of signature panels
const TOP_PANELS_PX = 300;
const SECTION_GAP_PX = 70;
const SIGNATURE_PANEL_PX = 60;
const MIN_SIGNATURES_PX = 240;
const LABEL_BOX_PX = 21;
const MARGIN_TOP = 150;
const MARGIN_BOTTOM = 90;

const baseYAxis = {
  linecolor: '#D3D3D3',
  linewidth: 1,
  ticks: '',
  mirror: 'all',
  tickfont: { family: 'Arial' },
};

const axisId = (i) => (i === 1 ? 'y' : `y${i}`);
const axisKey = (i) => (i === 1 ? 'yaxis' : `yaxis${i}`);
const center = ([y0, y1]) => (y0 + y1) / 2;

function groupStart(groups, index) {
  return groups
    .slice(0, index)
    .reduce((total, group) => total + group.data.length, 0);
}

function groupBorders(groups, y0, y1) {
  return groups.map((_, i) => ({
    type: 'rect',
    xref: 'x',
    yref: 'paper',
    x0: groupStart(groups, i) - 0.5,
    x1: groupStart(groups, i + 1) - 0.5,
    y0,
    y1,
    line: { width: 1 },
  }));
}

function groupBarTraces(groups, colors, axes) {
  return groups.map((group, i) => ({
    name: group.mutation,
    type: 'bar',
    marker: { color: colors[group.mutation] },
    x: group.data.map((_, j) => j + groupStart(groups, i)),
    y: group.data.map((e) => e.mutations || e.contribution || 0),
    hoverinfo: 'x+y',
    showlegend: false,
    ...axes,
  }));
}

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
  const seqmatrixData = data[2].data;

  // --- exposures ---
  const exposureBySignature = groupBy(
    exposureData?.filter((e) => e.exposure > 0.01),
    'signatureName'
  );
  const signatureNames = Object.keys(exposureBySignature);

  const longest = Math.max(0, ...signatureNames.map((e) => e.length));
  const extraMargin = longest < 7 ? 200 : longest * 12.5;

  const isSBS = arrayContainsTerms(signatureNames, ['SBS']);
  const signatureColors = isSBS ? colorPallet : colorPallet1;
  const signatureColor = (name, i) =>
    isSBS
      ? signatureColors[name.replace(/^\D*/, '').replace(')', '')]
      : signatureColors[i];

  const exposureSum = Object.values(exposureBySignature)
    .flat()
    .reduce((n, { exposure }) => n + exposure, 0);

  const percentSignature = Object.values(exposureBySignature).map(([e]) => ({
    signatureName: e.signatureName,
    percent: e.exposure / exposureSum,
  }));
  const percentByName = new Map(
    percentSignature.map((e) => [e.signatureName, e.percent])
  );

  const ptext = percentSignature
    .map((e) => `${(e.percent * 100).toFixed(1)}%*${e.signatureName}`)
    .join(' + ');

  // --- original profile, restricted to mutation types present in the signatures ---
  const signatureRows = Object.values(
    groupBy(
      signatureData?.filter((e) => percentByName.has(e.signatureName)),
      'signatureName'
    )
  );
  const signatureMutationTypes = new Set(
    signatureRows.flat().map((e) => e.mutationType)
  );
  const seqmatrixDataFilter = (seqmatrixData ?? []).filter((e) =>
    signatureMutationTypes.has(e.mutationType)
  );

  const mutationOrder = Object.keys(colors);
  const mutationGroupSort = (a, b) =>
    mutationOrder.indexOf(a.mutation) - mutationOrder.indexOf(b.mutation);

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

  const groupOriginal = groupFilterDataByMutation(
    normalizedOriginal,
    mutationRegex,
    mutationGroupSort
  );

  // signature/deconstructed rows aren't guaranteed to arrive in the same order as the original profile
  const mutationTypeIndex = new Map(
    groupOriginal
      .flatMap((g) => g.data.map((e) => e.mutationType))
      .map((m, i) => [m, i])
  );
  const byOriginalMutationTypeOrder = (a, b) =>
    mutationTypeIndex.get(a.mutationType) -
    mutationTypeIndex.get(b.mutationType);

  // --- deconstructed profile: exposure-weighted sum of signatures ---
  const deconstructed = Object.values(
    signatureRows.flat().reduce((acc, e) => {
      acc[e.mutationType] ??= { mutationType: e.mutationType, mutations: 0 };
      acc[e.mutationType].mutations +=
        e.contribution * percentByName.get(e.signatureName);
      return acc;
    }, {})
  );
  const deconstructedByType = new Map(
    deconstructed.map((e) => [e.mutationType, e.mutations])
  );

  const groupDeconstructed = groupDataByMutation(
    deconstructed,
    mutationRegex,
    mutationGroupSort,
    byOriginalMutationTypeOrder
  );

  const maxMutations = Math.max(
    getMaxMutations(normalizedOriginal) /
      getTotalMutations(normalizedOriginal),
    getMaxMutations(deconstructed) / getTotalMutations(deconstructed)
  );

  // --- per-signature profiles ---
  const groupSamples = signatureRows.map((rows) =>
    groupDataByMutation(
      rows,
      mutationRegex,
      mutationGroupSort,
      byOriginalMutationTypeOrder
    )
  );
  const signatureCount = groupSamples.length;

  // --- geometry (paper coordinates) ---
  const signaturesPx = Math.max(
    signatureCount * SIGNATURE_PANEL_PX,
    MIN_SIGNATURES_PX
  );
  const plotAreaPx = TOP_PANELS_PX + SECTION_GAP_PX + signaturesPx;
  const px = (pixels) => pixels / plotAreaPx;
  const gapSm = px(8);
  const gapLg = px(17);
  const labelBoxHeight = px(LABEL_BOX_PX);

  const signaturesTop = px(signaturesPx);
  const divide2 = signatureCount ? signaturesTop / signatureCount : 0;
  const divide1 = px(TOP_PANELS_PX) / 3;

  const signatureDomain = (i) => [
    i === 0 ? 0 : divide2 * i - gapSm,
    divide2 * (i + 1) - gapLg,
  ];
  // 0 = original, 1 = deconstructed, 2 = difference
  const topDomain = (k) => [
    1 - divide1 * (k + 1) - gapSm,
    k === 0 ? 1 : 1 - divide1 * k - gapLg,
  ];

  const topAxisBase = Math.max(signatureCount, 1);
  const [differenceAxis, deconstructedAxis, originalAxis] = [1, 2, 3].map(
    (k) => topAxisBase + k
  );

  // --- traces ---
  const originalTraces = groupBarTraces(groupOriginal, colors, {
    xaxis: 'x2',
    yaxis: axisId(originalAxis),
  });
  const deconstructedTraces = groupBarTraces(groupDeconstructed, colors, {
    xaxis: 'x2',
    yaxis: axisId(deconstructedAxis),
  });
  const differenceTraces = originalTraces.map((trace, i) => ({
    ...trace,
    y: groupOriginal[i].data.map(
      (e, j) => trace.y[j] - (deconstructedByType.get(e.mutationType) || 0)
    ),
    yaxis: axisId(differenceAxis),
  }));
  const signatureTraces = groupSamples.flatMap((groups, i) =>
    groupBarTraces(groups, colors, { yaxis: axisId(i + 1) })
  );

  const originalValues = originalTraces.flatMap((t) => t.y);
  const deconstructedValues = groupOriginal.flatMap((g) =>
    g.data.map((e) => deconstructedByType.get(e.mutationType) || 0)
  );
  const rss = getRss(differenceTraces.flatMap((t) => t.y));
  const cosineSimilarity = getCosineSimilarity(
    originalValues,
    deconstructedValues
  );

  const traces = [
    ...signatureTraces,
    ...differenceTraces,
    ...originalTraces,
    ...deconstructedTraces,
  ];

  const title = {
    text:
      '<b>Mutational Signature Association</b><br><b>RSS = ' +
      rss +
      '; Cosine Similarity = ' +
      cosineSimilarity +
      '</b>',
    font: { family: 'Times New Roman', size: 20 },
  };

  if (!traces.length) {
    return {
      traces,
      layout: {
        hoverlabel: { bgcolor: '#FFF' },
        height: 250,
        autosize: true,
        title,
        xaxis: {
          autorange: true,
          showgrid: false,
          zeroline: false,
          showline: false,
          ticks: '',
          showticklabels: false,
        },
        yaxis: {
          autorange: true,
          showgrid: false,
          zeroline: false,
          showline: false,
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
            font: { size: 28, color: 'grey' },
          },
        ],
      },
    };
  }

  // --- shapes ---
  const labelBoxes = (groups, y0) =>
    groupBorders(groups, y0, y0 + labelBoxHeight).map((shape, i) => ({
      ...shape,
      fillcolor: colors[groups[i].mutation],
    }));

  const sortedPercents = [...percentSignature].sort(
    (a, b) => b.percent - a.percent
  );
  let cumulative = 0;
  const percentBands = sortedPercents.map((e, i) => {
    const y0 = cumulative;
    cumulative += e.percent;
    return {
      name: e.signatureName,
      y0,
      y1: cumulative,
      color: signatureColor(e.signatureName, i),
    };
  });

  const percentShapes = percentBands.flatMap(({ y0, y1, color }) => [
    {
      type: 'rect',
      xref: 'paper',
      yref: 'paper',
      x0: -0.045,
      x1: -0.055,
      y0,
      y1,
      fillcolor: color,
      line: { width: 0 },
    },
    {
      type: 'line',
      xref: 'paper',
      yref: 'paper',
      x0: -0.055,
      x1: -0.058,
      y0: (y0 + y1) / 2,
      y1: (y0 + y1) / 2,
      line: { width: 1, color },
    },
  ]);

  const shapes = [
    ...(groupSamples[0] ? labelBoxes(groupSamples[0], signaturesTop - gapLg) : []),
    ...labelBoxes(groupOriginal, 1),
    ...groupSamples.flatMap((groups, i) =>
      groupBorders(groups, ...signatureDomain(i))
    ),
    ...[0, 1, 2].flatMap((k) => groupBorders(groupOriginal, ...topDomain(k))),
    ...percentShapes,
  ];

  // --- annotations ---
  const labelAnnotations = (groups, y) =>
    groups.map((group, i) => ({
      xref: 'x',
      yref: 'paper',
      xanchor: 'center',
      yanchor: 'bottom',
      x: groupStart(groups, i) + (group.data.length - 1) / 2,
      y,
      text: formatMutationLabels(group),
      showarrow: false,
      font: { color: 'white', size: 13, family: 'Times New Roman' },
      align: 'center',
    }));

  const sideLabel = (text, y) => ({
    xref: 'paper',
    yref: 'paper',
    xanchor: 'center',
    yanchor: 'middle',
    align: 'center',
    x: 1.017,
    y,
    text,
    textangle: 90,
    showarrow: false,
  });

  const axisTitle = (y) => ({
    xref: 'paper',
    yref: 'paper',
    xanchor: 'right',
    yanchor: 'middle',
    align: 'center',
    x: -0.025,
    y,
    text: '<b>Relative contribution</b>',
    font: { size: 16, family: 'Times New Roman' },
    textangle: -90,
    showarrow: false,
  });

  const annotations = [
    ...(groupSamples[0]
      ? labelAnnotations(groupSamples[0], signaturesTop - gapLg + px(2))
      : []),
    ...labelAnnotations(groupOriginal, 1 + px(2)),
    ...groupSamples.map((groups, i) =>
      sideLabel(groups[0].data[0].signatureName, center(signatureDomain(i)))
    ),
    sideLabel('Original', center(topDomain(0))),
    sideLabel('Deconstructed', center(topDomain(1))),
    sideLabel('Difference', center(topDomain(2))),
    axisTitle(center(topDomain(1))),
    axisTitle(signaturesTop / 2),
    {
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      yanchor: 'top',
      x: 0.5,
      y: 0,
      yshift: -45,
      text: '<b>Original Profile = ' + ptext + '</b>',
      font: { size: 13 },
      showarrow: false,
      align: 'center',
    },
    ...percentBands.map(({ name, y0, y1, color }) => ({
      xref: 'paper',
      yref: 'paper',
      xanchor: 'right',
      yanchor: 'middle',
      x: -0.058,
      y: (y0 + y1) / 2,
      text: name,
      font: { color },
      showarrow: false,
      align: 'center',
    })),
  ];

  // --- axes ---
  const tickLabels = groupOriginal.length ? formatTickLabels(groupOriginal) : [];
  const baseXAxis = {
    showline: true,
    tickangle: tickAngle,
    tickfont: { family: 'Courier New, monospace' },
    tickmode: 'array',
    tickvals: tickLabels.map((_, i) => i),
    ticktext: tickLabels,
    linecolor: '#E0E0E0',
    linewidth: 1,
    mirror: 'all',
    ticks: '',
  };

  const signatureAxes = Object.fromEntries(
    groupSamples.map((_, i) => [
      axisKey(i + 1),
      {
        ...baseYAxis,
        autorange: true,
        showgrid: true,
        gridcolor: '#F5F5F5',
        domain: signatureDomain(i),
        anchor: 'x',
      },
    ])
  );

  const differenceMax = findMaxAbsoluteYValue(differenceTraces);

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: plotAreaPx + MARGIN_TOP + MARGIN_BOTTOM,
    autosize: true,
    title,
    xaxis: { ...baseXAxis, anchor: 'y' },
    xaxis2: { ...baseXAxis, anchor: axisId(differenceAxis) },
    ...signatureAxes,
    [axisKey(differenceAxis)]: {
      ...baseYAxis,
      autorange: false,
      range: [-differenceMax * 1.5, differenceMax * 1.5],
      showgrid: true,
      gridcolor: '#F5F5F5',
      domain: topDomain(2),
      anchor: 'x2',
    },
    [axisKey(deconstructedAxis)]: {
      ...baseYAxis,
      autorange: false,
      range: [0, maxMutations * 1.2],
      domain: topDomain(1),
      anchor: 'x2',
    },
    [axisKey(originalAxis)]: {
      ...baseYAxis,
      autorange: false,
      range: [0, maxMutations * 1.2],
      domain: topDomain(0),
      anchor: 'x2',
    },
    shapes,
    annotations,
    margin: { l: extraMargin, t: MARGIN_TOP, b: MARGIN_BOTTOM },
  };

  return { traces, layout };
}
