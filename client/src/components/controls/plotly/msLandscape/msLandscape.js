import { groupBy } from 'lodash';

export default function MsLandscape(data, arg) {
  console.log(data);
  console.log(arg);
  const rawDataActivity = data[0]['data']; //exposure data
  const rawDataSignature = data[1]['data']; //signature data
  const rawDataSpectrum = data[2]['data']; //seqmatrix data
  console.log('activity/exposure data');
  console.log(rawDataActivity);
  console.log('signature data');
  console.log(rawDataSignature);
  console.log('spectrum/seqmatrix data');
  console.log(rawDataSpectrum);

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
  // const heatmapColorscale = [
  //   [0, 'rgb(34,7,139)'],
  //   [0.1, 'rgb(34,7,139)'],
  //   [0.2, 'rgb(34,7,139)'],
  //   [0.3, 'rgb(34,7,139)'],
  //   [0.4, 'rgb(34,7,139)'],
  //   [0.5, 'rgb(34,7,139)'],
  //   [0.6, 'rgb(34,7,139)'],
  //   [0.7, 'rgb(34,7,139)'],
  //   [0.7, 'rgb(128,12,162)'],
  //   [0.8, 'rgb(200,69,125)'],
  //   [0.8, 'rgb(208,79,117)'],
  //   [0.9, 'rgb(245,149,70)'],
  //   [0.9, 'rgb(247,152,67)'],
  //   [1, 'rgb(241,246,34)'],
  // ];
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

    [0.6, 'rgb(128,12,162)'],
    [0.7, 'rgb(128,12,162)'],

    [0.7, 'rgb(200,69,125)'],
    [0.8, 'rgb(200,69,125)'],

    [0.8, 'rgb(208,79,117)'],
    [0.9, 'rgb(245,149,70)'],

    [0.9, 'rgb(247,152,67)'],
    [1.0, 'rgb(241,246,34)'],
  ];

  const groupBySample_activity = groupBy(rawDataActivity, 'sample');
  console.log(groupBySample_activity);

  const groupBySignatureName_activity = groupBy(
    rawDataActivity,
    'signatureName'
  );
  console.log(groupBySignatureName_activity);

  // const ordered = {};
  // Object.keys(groupBySignatureName_activity)
  //   .sort()
  //   .reverse()
  //   .forEach(function (key) {
  //     ordered[key] = groupBySignatureName_activity[key];
  //   });

  // console.log(ordered);

  const dataSignature = Object.entries(groupBySample_activity).map(
    ([key, value]) => ({
      signatureName: key,
      total: value.reduce((a, e) => a + parseInt(e.exposure), 0),
      sample: value.map((e) => e.sample),
      exposure: value.map((e) => e.exposure),
    })
  );
  console.log(dataSignature);
  console.log(dataSignature[1]);

  const groupByMutationType_signature = groupBy(
    rawDataSignature,
    'mutationType'
  );
  console.log(groupByMutationType_signature);

  const groupByMutationType_spectrum = groupBy(rawDataSpectrum, 'mutationType');
  console.log(groupByMutationType_spectrum);
  console.log(groupBy(rawDataSpectrum, 'cancer'));
  const orignal_genomes = rawDataSpectrum;
  //const signature = groupByMutationType_signature;
  const signature = rawDataSignature;
  const signature_activaties = rawDataActivity;

  //  decompsite_input <- calculate_similarities(orignal_genomes = seqmatrix_refdata_input, signature = signature_refsets_input, signature_activaties = exposure_refdata_input)
  //calculate_similarities(orignal_genomes = seqmatrix_refdata_input, signature = signature_refsets_input, signature_activaties = exposure_refdata_input)
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
  //const cosine = cosineSimilarity(s1mutations, s2mutations).toFixed(3);

  function uniq(a) {
    return Array.from(new Set(a));
  }
  function calculate_similarities(
    orignal_genomes,
    signature,
    signature_activaties
  ) {
    const data2 = orignal_genomes;
    const data3 = signature;
    const data4 = signature_activaties;
    console.log(data2);
    console.log(data3);
    console.log(data4);

    const sample_name_total = data4.map((e) => e.sample);
    const sample_name_tmp = data2.map((e) => e.sample);
    const sample_name = uniq(
      sample_name_total.filter((x) => sample_name_tmp.includes(x))
    );
    console.log(sample_name);

    const signature_name_total = data4.map((e) => e.signatureName);
    const signature_name_tmp = data3.map((e) => e.signatureName);
    const signature_name = uniq(
      signature_name_total.filter((x) => signature_name_tmp.includes(x))
    );
    console.log(signature_name_total);
    console.log(signature_name_tmp);
    console.log(signature_name);
  }

  calculate_similarities(orignal_genomes, signature, signature_activaties);

  const barCharColor = Object.keys(groupBySignatureName_activity);
  console.log(barCharColor);

  const xAxisName = Object.values(groupBySignatureName_activity)[0]
    .map((e) => ({
      sample: e.sample,
    }))
    .flat();
  console.log(xAxisName.map((e) => e));
  const traces1 = Object.entries(groupBySignatureName_activity)
    .reverse()
    .map(([key, value]) => ({
      key: key,
      value: value,
      map: barCharColor.map((e) => e),
      co: barCharColor.map((e) => e.replace(/^\D*/, '')),
      //showlegend: true,
      name: key,
      type: 'bar',
      x: value.map((e) => e.sample),
      y: value.map((e, i) => e.exposure / dataSignature[i].total),
      total: value.map((e, i) => dataSignature[i].total),
      marker: {
        color: colors[key.replace(/^\D*/, '')],
      },
      showlegend: false,
      test: value.map((d, i) => d.exposure / dataSignature[i].total),
    }));
  console.log(traces1);
  const traces2 = [
    {
      z: [xAxisName.map((_, i) => i / 65)],
      x: xAxisName.map((e) => e.sample),
      hoverongaps: false,
      xaxis: 'x',
      yaxis: 'y2',
      type: 'heatmap',
      colorscale: heatmapColorscale,
      colorbar: {
        title: { text: 'cosine similarity' },
        orientation: 'h',
        //borderwidth: 1,
        //xanchor: 'left',
        //yanchor: 'top',
        //y: 1,
        tickmode: 'array',
        tickvals: [0, 0.6, 0.7, 0.8, 0.9, 1],
        //len: 0.15,
        //thickness: 10,
        //x: 0,
      },
      xgap: 0.2,
    },
  ];
  console.log(traces2);

  const traces3 = Object.entries(groupBySignatureName_activity)
    .reverse()
    .map(([key, value]) => ({
      key: key,
      value: value,
      map: barCharColor.map((e) => e),
      co: barCharColor.map((e) => e.replace(/^\D*/, '')),
      name: key,
      type: 'bar',
      x: value.map((e) => e.sample),
      y: value.map((e, i) => e.exposure),
      //showlegend: true,
      marker: {
        color: colors[key.replace(/^\D*/, '')],
      },
      xaxis: 'x',
      yaxis: 'y3',
      // transforms: [
      //   {
      //     type: 'sort',
      //     target: 'y',
      //     order: 'descending',
      //   },
      // ],
    }));
  console.log(traces3);

  const traces = [...traces3, ...traces2, ...traces1];
  console.log(traces);

  const shapes = [];

  const lines = [];

  let annotations = [];

  const layout = {
    autosize: true,
    height: 1080,
    barmode: 'stack',
    hovermode: 'closest',
    legend: {
      orientation: 'h',
      // title: { text: 'Signatures Name' },
      traceorder: 'reversed',
    },

    xaxis: {
      tickmode: 'array',
      tickvals: xAxisName.map((_, i) => i),
      ticktext: xAxisName.map((e) => e.sample),
      type: 'category',
      tickangle: -90,
    },
    yaxis: { title: 'Signature contribution', domain: [0, 0.45] },
    yaxis2: { title: '', domain: [0.46, 0.52] },
    yaxis3: { title: 'Number of mutation', domain: [0.55, 1] },

    shapes: [...shapes, ...lines],
    annotations: annotations,
  };

  var config = {
    //responsive: true,
  };

  return { traces, layout, config };
}
