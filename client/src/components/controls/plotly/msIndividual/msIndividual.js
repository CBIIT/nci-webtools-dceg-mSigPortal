import { groupBy } from 'lodash';
import {
  groupDataByMutation,
  getTotalMutations,
  getMaxMutations,
  getCosineSimilarity,
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
    exposureData.filter((o) => o['exposure'] > 0),
    'signatureName'
  );
  console.log(exposure_groupBySignature);

  const signatureNames = Object.keys(exposure_groupBySignature).map((e) => e);
  console.log(signatureNames);
  const exposureSum = Object.values(exposure_groupBySignature)
    .flat()
    .reduce((n, { exposure }) => n + exposure, 0);
  console.log(exposureSum);

  const signature_groupBySignature = groupBy(
    signatureData.filter((e) => signatureNames.includes(e.signatureName)),
    'signatureName'
  );
  console.log(signature_groupBySignature);

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
  console.log(mutationTypes);

  const seqmatrix_groupByMutationType = groupBy(
    segmatrixData.filter((e) => mutationTypes.includes(e.mutationType)),
    'mutationType'
  );
  console.log(seqmatrix_groupByMutationType);

  const seqmatrixDataFilter = Object.values(
    seqmatrix_groupByMutationType
  ).flat(); //original data for the comparison
  console.log(seqmatrixDataFilter);

  const arraySignatureData = Object.values(signature_groupBySignature).map(
    (e) => e
  );

  console.log(arraySignatureData);
  let totalMutationsArray = [];
  let maxMutationsArray = [];
  for (let i = 0; i < arraySignatureData.length; i++) {
    console.log(arraySignatureData[i]);
    totalMutationsArray.push(getTotalMutations(arraySignatureData[i]));
    maxMutationsArray.push(
      getMaxMutations(arraySignatureData[i]) /
        getTotalMutations(arraySignatureData[i])
    );
  }
  const maxMutations = Math.max(...maxMutationsArray);
  console.log(totalMutationsArray);
  console.log(maxMutationsArray);
  console.log(maxMutations);

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
  const traces = [];
  for (var i = 0; i < groupSamples.length; i++) {
    for (var j = 0; j < groupSamples[i].length; j++) {
      let t = {
        name: groupSamples[i][j].mutation,
        type: 'bar',
        maker: { color: colors[groupSamples[i][j].mutation] },
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
      traces.push(t);
    }
  }
  const mutationAnnotation = groupSamples[0].map(
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
      y: 1.005,
      text: formatMutationLabels(group),
      showarrow: false,
      font: { color: 'white' },
      align: 'center',
    })
  );

  const mutationLabelBox = groupSamples[0].map((group, groupIndex, array) => ({
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
  console.log(traces);
  const sampleLabels = [];

  const tickLabels = formatTickLabels(groupSamples[0]);
  console.log(tickLabels);

  const layout = {
    showlegend: true,
    hoverlabel: { bgcolor: '#FFF' },
    height: 900,

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
      //range: [0, maxMutations * 1.2],
      linecolor: '#D3D3D3',
      linewidth: 1,
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      showgrid: true,
      gridcolor: '#F5F5F5',
      domain: [0, 0.18],
    },
    yaxis2: {
      autorange: true,
      // range: [0, maxMutations * 1.2],
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      domain: [0.2, 0.38],
      title: { text: '<b>Relative contribution</b>' },
    },
    yaxis3: {
      autorange: true,
      // range: [0, maxMutations * 1.2],
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      domain: [0.4, 0.58],
    },
    yaxis4: {
      autorange: true,
      // range: [0, maxMutations * 1.2],
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      domain: [0.6, 0.78],
    },
    yaxis5: {
      autorange: true,
      // range: [0, maxMutations * 1.2],
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      tickfont: {
        family: 'Arial',
      },
      domain: [0.8, 0.88],
    },
    shapes: [...mutationLabelBox],
    annotations: [
      ...mutationAnnotation,
      sampleLabels,
      //yTitleAnnotation,
    ],
  };

  console.log(layout);
  return { traces: traces, layout: layout };
}
