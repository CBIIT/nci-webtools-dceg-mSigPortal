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

  const groupBySample_activity = groupBy(rawDataActivity, 'sample');
  console.log(groupBySample_activity);

  const groupBySignatureName_activity = groupBy(
    rawDataActivity,
    'signatureName'
  );
  console.log(groupBySignatureName_activity);

  const dataSignature = Object.entries(groupBySignatureName_activity).map(
    ([key, value]) => ({
      signatureName: key,
      data: value,
    })
  );
  console.log(dataSignature);

  const groupByMutationType_signature = groupBy(
    rawDataSignature,
    'mutationType'
  );
  console.log(groupByMutationType_signature);

  const groupByMutationType_spectrum = groupBy(rawDataSpectrum, 'mutationType');
  console.log(groupByMutationType_spectrum);
  console.log(groupBy(rawDataSpectrum, 'cancer'));
  const orignal_genomes = rawDataSpectrum;
  const signature = groupByMutationType_signature;
  const signature_activaties = rawDataActivity;

  //  decompsite_input <- calculate_similarities(orignal_genomes = seqmatrix_refdata_input, signature = signature_refsets_input, signature_activaties = exposure_refdata_input)
  //calculate_similarities(orignal_genomes = seqmatrix_refdata_input, signature = signature_refsets_input, signature_activaties = exposure_refdata_input)

  function calculate_similarities(
    orignal_genomes,
    signature,
    signature_activaties
  ) {
    const data2 = orignal_genomes;
    const data3 = signature;
    const data4 = signature_activaties;
    console.log(data2);
  }
  const barCharColor = Object.keys(groupBySignatureName_activity);
  console.log(barCharColor);

  const traces = Object.entries(groupBySignatureName_activity).map(
    ([key, value]) => ({
      key: key,
      value: value,
      map: barCharColor.map((e) => e),
      co: barCharColor.map((e) => e.replace(/^\D*/, '')),
      name: key,
      type: 'bar',
      x: value.map((e) => e.sample),
      y: value.map((e) => e.exposure),
      // marker: {
      //   color: barCharColor.map((e) => colors[e.replace(/^\D*/, '')]),
      // },
    })
  );

  // const traces = Object.entries(groupBySample_activity).map(([key, value]) => ({
  //   key: key,
  //   value: value,
  //   map: value.map((group) => group.signatureName),
  //   name: key,
  //   type: 'bar',
  //   x: xValue,
  //   y: value.map((e) => e.exposure),
  //   hovertemplate: '%{y}, %{x} <extra></extra>',
  //   marker: {
  //     colors: value.map(
  //       (group) =>
  //         colors[group.signatureName.replace(/^\D*/, '').replace(')', '')]
  //     ),
  //   },
  // }));
  console.log(traces);

  const shapes = [];

  const lines = [];

  let annotations = [];

  const layout = {
    // width: totalCancer > 1 ? null : 350,
    autosize: true,
    height: 500,
    barmode: 'stack',

    yaxis: { title: 'Number of mutation' },

    shapes: [...shapes, ...lines],
    annotations: annotations,
  };

  var config = {
    //responsive: true,
  };

  return { traces, layout, config };
}
