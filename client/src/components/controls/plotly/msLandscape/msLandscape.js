import { groupBy } from 'lodash';
import { First } from 'react-bootstrap/esm/PageItem';

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

  // const ordered = {};
  // Object.keys(groupBySignatureName_activity)
  //   .sort()
  //   .reverse()
  //   .forEach(function (key) {
  //     ordered[key] = groupBySignatureName_activity[key];
  //   });

  // console.log(ordered);
  const mutationType = [];

  const mType = groupBy(rawDataSignature, 'mutationType');
  Object.keys(mType).map((e) => mutationType.push(e));
  console.log(mutationType);

  const filterDataSpectrum = rawDataSpectrum.filter((e) =>
    mutationType.includes(e.mutationType)
  );

  const filterDataSignature = rawDataSignature.filter((e) =>
    mutationType.includes(e.mutationType)
  );
  console.log(filterDataSpectrum);
  const exposure_refdata_input = rawDataActivity.map((group) => ({
    exposure: group.exposure,
    sample: group.sample,
    signatureName: group.signatureName,
  }));
  //.filter((obj) => obj.exposure !== 0);
  console.log('exposure_refdata_input');
  console.log(exposure_refdata_input);

  const seqmatrix_refdata_input = filterDataSpectrum.map((group) => ({
    mutation: group.mutations,
    mutationType: group.mutationType,
    sample: group.sample,
  }));
  console.log('seqmatrix_refdata_input');
  console.log(seqmatrix_refdata_input);

  const signature_refsets_input = filterDataSignature.map((group) => ({
    contribution: group.contribution,
    mutationType: group.mutationType,
    signatureName: group.signatureName,
  }));

  console.log('signature_refsets_input');
  console.log(signature_refsets_input);

  const groupBySample_exposure = groupBy(exposure_refdata_input, 'sample');
  console.log(groupBySample_exposure);

  const groupBySignatureName_exposure = groupBy(
    exposure_refdata_input,
    'signatureName'
  );
  console.log(groupBySignatureName_exposure);

  function groupBy1(objectArray, property) {
    return objectArray.reduce(function (acc, obj) {
      var key = obj[property];
      if (!acc[key]) {
        acc[key] = [];
      }
      if (obj.exposure !== 0) {
        acc[key].push(obj);
      }

      return acc;
    }, {});
  }
  const groupBySignatureName_exposure2 = groupBy1(
    exposure_refdata_input,
    'signatureName'
  );
  console.log(groupBySignatureName_exposure2);

  const dataSignature = Object.entries(groupBySample_exposure).map(
    ([key, value]) => ({
      signatureName: key,
      total: value.reduce((a, e) => a + parseInt(e.exposure), 0),
      sample: value.map((e) => e.sample),
      exposure: value.map((e) => e.exposure),
    })
  );
  //.filter((obj) => obj.exposure !== 0);
  console.log(dataSignature);

  const groupByMutationType_signature = groupBy(
    rawDataSignature,
    'mutationType'
  );
  console.log(groupByMutationType_signature);

  const groupByMutationType_spectrum = groupBy(
    filterDataSpectrum,
    'mutationType'
  );
  console.log(groupByMutationType_spectrum);
  console.log(groupBy(filterDataSpectrum, 'cancer'));

  const sortSignatureName = (sourceArray) => {
    const sortByLocation = (a, b) =>
      a.signatureName.localeCompare(b.signatureName, 'en', { numeric: true });
    return sourceArray.sort(sortByLocation);
  };
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

  //matrix multiplication
  function matrixDot(A, B) {
    var result = new Array(A.length)
      .fill(0)
      .map((row) => new Array(B[0].length).fill(0));

    return result.map((row, i) => {
      return row.map((val, j) => {
        return A[i].reduce((sum, elm, k) => sum + elm * B[k][j], 0);
      });
    });
  }

  const multiplyMatrices = (a, b) => {
    if (!Array.isArray(a) || !Array.isArray(b) || !a.length || !b.length) {
      throw new Error('arguments should be in 2-dimensional array format');
    }
    let x = a.length,
      z = a[0].length,
      y = b[0].length;
    console.log(x);
    console.log(y);
    console.log(z);
    console.log(b.length);
    if (b.length !== z) {
      // XxZ & ZxY => XxY
      throw new Error(
        'number of columns in the first matrix should bethe same as the number of rows in the second'
      );
    }
    let productRow = Array.apply(null, new Array(y)).map(
      Number.prototype.valueOf,
      0
    );
    let product = new Array(x);
    for (let p = 0; p < x; p++) {
      product[p] = productRow.slice();
    }
    for (let i = 0; i < x; i++) {
      for (let j = 0; j < y; j++) {
        for (let k = 0; k < z; k++) {
          product[i][j] += a[i][k] * b[k][j];
        }
      }
    }
    return product;
  };

  function matrixMul(m1, m2) {
    const fil_m1 = m1.length;
    const col_m1 = m1[0].length;
    const fil_m2 = m2.length;
    const col_m2 = m2[0].length;
    console.log('m1 lenght: ' + fil_m1);
    console.log('m1 col: ' + col_m1);
    console.log('m2 lenght: ' + fil_m2);
    console.log('m2 col: ' + col_m2);
    if (col_m1 != fil_m2) throw 'Matrices cannot be multiplied';
    let multiplication = new Array(fil_m1);
    for (let x = 0; x < multiplication.length; x++)
      multiplication[x] = new Array(col_m2).fill(0);
    for (let x = 0; x < multiplication.length; x++) {
      for (let y = 0; y < multiplication[x].length; y++) {
        for (let z = 0; z < col_m1; z++) {
          multiplication[x][y] = multiplication[x][y] + m1[x][z] * m2[z][y];
        }
      }
    }
  }

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
    console.log('data2: ');
    console.log(data2);
    console.log('data3: ');
    console.log(data3);
    console.log('data4: ');
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
    //console.log(signature_name_total);
    //console.log(signature_name_tmp);
    console.log(signature_name);
    const newData3 = data3.filter((e) =>
      signature_name.includes(e.signatureName)
    );
    console.log(newData3);
    const dataSample2 = groupBy(
      data2.map((e) => ({
        mutation: e.mutation,
        sample: e.sample,
      })),
      'sample'
    );
    const dataSample3 = groupBy(
      sortSignatureName(newData3).map((e) => ({
        contribution: e.contribution,
        signatureName: e.signatureName,
      })),
      'signatureName'
    );

    const dataSample4 = groupBy(
      sortSignatureName(data4).map((e) => ({
        exposure: e.exposure,
        signatureName: e.signatureName,
      })),
      'signatureName'
    );
    const array2 = [];
    Object.values(dataSample2).forEach((group) => {
      //console.log(group);
      array2.push(group.map((e) => e.mutation));
    });
    const array3 = [];
    Object.values(dataSample3).forEach((group) => {
      array3.push(group.map((e) => e.contribution));
    });
    const array4 = [];
    Object.values(dataSample4).forEach((group) => {
      array4.push(group.map((e) => e.exposure));
    });
    console.log(dataSample2);
    console.log(array2);
    console.log(dataSample3);
    console.log(array3);
    console.log(dataSample4);
    console.log(array4);
    const array3_trans = array3[0].map((_, colIndex) =>
      array3.map((row) => row[colIndex])
    );
    console.log(array3_trans);
    const genomes = array2;
    console.log(genomes);
    const est_genomes = multiplyMatrices(array3_trans, array4);
    console.log(est_genomes);
  }

  calculate_similarities(
    seqmatrix_refdata_input,
    signature_refsets_input,
    exposure_refdata_input
  );

  const barCharColor = Object.keys(groupBySignatureName_exposure);
  console.log(barCharColor);

  const barCharColor2 = Object.keys(groupBySignatureName_exposure2);
  console.log(barCharColor2);

  const xAxisName = Object.values(groupBySignatureName_exposure)[0]
    .map((e) => ({
      sample: e.sample,
    }))
    .flat();
  console.log(xAxisName.map((e) => e));
  const traces1 = Object.entries(groupBySignatureName_exposure)
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
      expo: value.map((e, i) => e.exposure),
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
      zmin: 0.6,
      zmax: 1,
      colorbar: {
        orientation: 'h',
        x: 0.5,
        y: -0.3,
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
      xgap: 0.2,
    },
  ];
  console.log(traces2);

  const traces3 = Object.entries(groupBySignatureName_exposure)
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

  const longest = xAxisName.reduce(
    (a, e) => (a > e.sample.length ? a : e.sample.length),
    0
  );
  const extraMargin = longest < 10 ? -0.157 : (longest * -0.027) / 2;
  console.log(longest);
  console.log(extraMargin);

  const text = {
    x: 0,
    //y: -0.157,
    y: extraMargin,
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
      // title: { text: 'Signatures Name' },
      traceorder: 'reversed',
      x: 0,
      y: extraMargin,
    },

    xaxis: {
      tickmode: 'array',
      tickvals: xAxisName.map((_, i) => i),
      ticktext: xAxisName.map((e) => e.sample),
      type: 'category',
      tickangle: -90,
      ticks: '',
    },
    yaxis: { title: 'Signature contribution', domain: [0, 0.47] },
    yaxis2: {
      title: '',
      domain: [0.47, 0.51],
      showticklabels: false,
      ticks: '',
    },
    yaxis3: { title: 'Number of mutation', domain: [0.53, 1] },

    bargap: 0.04,
    shapes: [...shapes, ...lines],
    annotations: annotations,
  };

  var config = {
    //responsive: true,
  };

  return { traces, layout, config };
}
