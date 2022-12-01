import { groupBy } from 'lodash';
import { First } from 'react-bootstrap/esm/PageItem';

export default function MsLandscape(cosineData, exposureData) {
  console.log(cosineData);
  console.log(exposureData);

  let colors = {
    1: '#4a9855',
    2: '#e2a8ab',
    3: '#40004b',
    4: '#5aa1ca',
    5: '#305d39',
    6: '#785940',
    '7a': '#6e70b7',
    '7b': '#ff7f00',
    '7c': '#fec44f',
    '7d': '#846a2a',
    8: '#cab2d6',
    9: '#f4a582',
    '10a': '#8dd3c7',
    '10b': '#5e4fa2',
    '10c': '#761429',
    11: '#9e0142',
    12: '#ffed6f',
    13: '#e41a1c',
    14: '#ffffbf',
    15: '#4d4d4d',
    16: '#513276',
    '17a': '#df4c7d',
    '17b': '#08519c',
    18: '#b3de69',
    19: '#dfc27d',
    20: '#b2182b',
    21: '#9ecae1',
    22: '#01665e',
    23: '#d53e4f',
    24: '#1c9099',
    25: '#35978f',
    26: '#ec7014',
    27: '#f46d43',
    28: '#de77ae',
    29: '#fdae61',
    30: '#d9d9d9',
    31: '#f781bf',
    32: '#dd1c77',
    33: '#b25d7e',
    34: '#fee08b',
    35: '#fc8d59',
    36: 'yellow',
    37: '#e6f598',
    38: '#abdda4',
    39: '#636363',
    40: '#b15928',
    41: '#fccde5',
    42: '#ae017e',
    43: '#66c2a5',
    44: '#8c6bb1',
    45: '#3288bd',
    46: '#e6f598',
    47: '#bababa',
    48: '#5e4fa2',
    49: '#40004b',
    50: '#762a83',
    51: '#9970ab',
    52: '#c2a5cf',
    53: '#e7d4e8',
    54: '#fcc5c0',
    55: '#d9f0d3',
    56: '#8c510a',
    57: '#a6dba0',
    58: '#5aae61',
    59: '#1b7837',
    60: '#00441b',
    84: '#063C3C',
    85: '#AA9139',
    92: '#0E1844',
    '-others': '#cececa',
  };

  const heatmapColorscale = [
    [0, 'rgb(34,7,139)'],
    [0.1, 'rgb(34,7,139)'],

    [0.1, 'rgb(34,7,139)'],
    [0.2, 'rgb(34,7,139)'],

    [0.2, 'rgb(34,7,139)'],
    [0.3, 'rgb(34,7,139)'],

    [0.3, 'rgb(34,7,139)'],
    [0.4, 'rgb(34,7,139)'],

    [0.4, 'rgb(34,7,139)'],
    [0.5, 'rgb(34,7,139)'],

    [0.5, 'rgb(34,7,139)'],
    [0.6, 'rgb(34,7,139)'],

    [0.6, 'rgb(34,7,139)'],
    [0.7, 'rgb(124,9,163)'],

    [0.7, 'rgb(145,22,156)'],
    [0.8, 'rgb(202,72,122)'],

    [0.8, 'rgb(220,94,103)'],
    [0.9, 'rgb(244,145,72)'],

    [0.9, 'rgb(252,165,55)'],
    [1.0, 'rgb(241,246,34)'],
  ];

  const sortSignatureName = (sourceArray) => {
    const sortByLocation = (a, b) =>
      a.signatureName.localeCompare(b.signatureName, 'en', { numeric: true });
    return sourceArray.sort(sortByLocation);
  };

  // const groupBySignatureName_exposure = groupBy(
  //   sortSignatureName(exposureData),
  //   'signatureName'
  // );
  // console.log(groupBySignatureName_exposure);

  const groupBySample_exposure = groupBy(exposureData, 'sample');
  console.log(groupBySample_exposure);

  const xAxis = Object.keys(groupBySample_exposure).flat();

  console.log(xAxis);

  const sortedCosin = cosineData.sort(
    (a, b) => xAxis.indexOf(a.sample) - xAxis.indexOf(b.sample)
  );
  console.log(sortedCosin);

  // const dataSignature = Object.entries(groupBySample_exposure).map(
  //   ([key, value]) => ({
  //     sample: key,
  //     total: value.reduce((a, e) => a + parseInt(e.exposure), 0),
  //     signatureName: value.map((e) => e.signatureName),
  //     exposure: value.map((e) => e.exposure),
  //     exposureNorm: value.map(
  //       (e) => e.exposure / value.reduce((a, e) => a + parseInt(e.exposure), 0)
  //     ),
  //   })
  // );
  // console.log(dataSignature);

  const dataSignature2 = Object.entries(groupBySample_exposure)

    .map(([sample, data]) => {
      const total = data.reduce((a, e) => a + parseInt(e.exposure), 0);

      return data.map((e) => ({ ...e, total }));
    })

    .flat();

  console.log(dataSignature2);

  const groupBySignatureName_exposure2 = groupBy(
    sortSignatureName(dataSignature2),
    //dataSignature2,
    'signatureName'
  );
  console.log(groupBySignatureName_exposure2);
  console.log(Object.entries(groupBySignatureName_exposure2).reverse());

  const barCharColor = Object.keys(groupBySignatureName_exposure2);
  console.log(barCharColor);

  const tracesNormalize = Object.entries(groupBySignatureName_exposure2)
    .reverse()
    .map(([key, value]) => ({
      key: key,
      value: value,
      map: barCharColor.map((e) => e),
      co: barCharColor.map((e) => e.replace(/^\D*/, '')),
      name: key,
      type: 'bar',
      x: value.map((e) => e.sample),
      y: value.map((e, i) => e.exposure / e.total),
      marker: {
        color: colors[key.replace(/^\D*/, '')],
      },
      showlegend: false,
      exposure: value.map((e, i) => e.exposure),
    }));
  console.log(tracesNormalize);
  const tracesHeatMap = [
    {
      z: [cosineData.map((e) => e.similarity)],
      x: cosineData.map((e) => e.sample),
      hoverongaps: false,
      xaxis: 'x',
      yaxis: 'y2',
      type: 'heatmap',
      colorscale: heatmapColorscale,
      zmin: 0.6,
      zmax: 1,
      colorbar: {
        orientation: 'h',
        x: 0.5,
        y: xAxis.length > 250 ? -0.2 : -0.3,
        bordercolor: 'gray',
        tickmode: 'array',
        tickvals: [0.6, 0.7, 0.8, 0.9, 1],
        title: {
          text: 'Cosine Similarity',
          font: {
            family: 'Arial',
            size: 17,
            color: 'rgb(37,37,37)',
          },
        },
      },
      xgap: 2,
    },
  ];
  console.log(tracesHeatMap);
  const tracesStackedBar = Object.entries(groupBySignatureName_exposure2)
    .reverse()
    .map(([key, value]) => ({
      key: key,
      value: value,
      map: barCharColor.map((e) => e),
      co: barCharColor.map((e) => e.replace(/^\D*/, '')),
      name: key,
      type: 'bar',
      x: value.filter((obj) => obj.exposure !== 0).map((e) => e.sample),
      y: value.filter((obj) => obj.exposure !== 0).map((e, i) => e.exposure),
      //showlegend: true,
      marker: {
        color: colors[key.replace(/^\D*/, '')],
      },
      xaxis: 'x',
      yaxis: 'y3',
      transforms: [
        {
          type: 'sort',
          target: 'y',
          order: 'descending',
        },
      ],
    }));
  console.log(tracesStackedBar);

  const traces = [...tracesHeatMap, ...tracesStackedBar, ...tracesNormalize];
  console.log(traces);

  const longest = xAxis.reduce((a, e) => (a > e.length ? a : e.length), 0);
  const extraMargin =
    longest > 0 && longest < 10 ? -0.157 : (longest * -0.027) / 2;
  console.log(longest);
  console.log(extraMargin);

  const text = {
    x: 0,
    y: xAxis.length > 250 ? -0.03 : extraMargin,

    xanchor: 'left',
    yanchor: 'bottom',
    xref: 'paper',
    yref: 'paper',
    text: 'Mutational Signatures:',
    showarrow: false,
    font: {
      family: 'Arial',
      size: 17,
      color: 'rgb(37,37,37)',
    },
  };
  console.log(text);
  const shapes = [];

  const lines = [];

  let annotations = [text];

  const layout = {
    autosize: true,
    height: 1200,
    barmode: 'stack',
    hovermode: 'closest',
    legend: {
      orientation: 'h',
      //title: { text: 'Signatures Name<br>' },
      traceorder: 'normal',
      x: 0,
      y: xAxis.length > 250 ? -0.03 : extraMargin,
    },

    xaxis: {
      tickmode: 'array',
      tickvals: xAxis.map((_, i) => i),
      ticktext: xAxis.map((e) => e),
      type: 'category',
      tickangle: -90,
      ticks: '',
      showticklabels: xAxis.length > 250 ? false : true,
      zeroline: false,
    },
    yaxis: { title: 'Signature contribution', domain: [0, 0.49] },
    yaxis2: {
      title: '',
      domain: [0.475, 0.493],
      showticklabels: false,
      ticks: '',
    },
    yaxis3: { title: 'Number of mutation', domain: [0.5, 1] },

    bargap: 0.04,
    shapes: [...shapes, ...lines],
    annotations: annotations,
  };

  var config = {
    //responsive: true,
  };

  return { traces, layout, config };
}
