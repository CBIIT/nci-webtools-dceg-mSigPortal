import { groupBy } from 'lodash';
import { First } from 'react-bootstrap/esm/PageItem';
import { get } from 'react-hook-form';
import CosineSimilarityPlot from '../../../pages/catalog/referenceSignature/cosineSimilarity/form-plot';
import { getRandomColor, isNumber, colorPallet0 } from '../../utils/utils';

export default function MsLandscape(cosineData, exposureData, variableData) {
  console.log(cosineData);
  console.log(exposureData);
  console.log(variableData);

  let charColors = {};
  let arrNames = [];
  let arrColors = [];
  let mixMatch = false;
  if (variableData.length > 0) {
    let stringData;
    if (!isNumber(variableData[0].value1)) {
      const stringVal = variableData.map((e) => e.value1);
      stringData = [...new Set(stringVal)].sort();
      //var bg_colors = stringData.map((e) => getRandomColor());
      var bg_colors = stringData.map(
        (e, i) => colorPallet0[i % colorPallet0.length]
      );
      console.log(stringData);
      console.log(bg_colors);
      stringData.forEach((element, index) => {
        charColors[element] = bg_colors[index];
      });
      arrNames = stringData;
      arrColors = bg_colors;
      if (
        variableData[0].hasOwnProperty('value2') &&
        isNumber(variableData[0].value2)
      ) {
        mixMatch = true;
      } else {
        mixMatch = false;
      }
    }
  }

  console.log(arrNames);
  console.log(arrColors);
  console.log(mixMatch);

  // const charColors = {
  //   A: '#4DBBD5',
  //   B: '#3C5488',
  //   C: '#E64B35',
  //   D: '#F39B7F',
  //   E: '#00A087',
  //   F: '#EBC6C4',
  // };

  const defaultNames = ['SBS', 'DBS', 'ID'];
  const names = exposureData.map((group) => group.signatureName);
  const contains = defaultNames.some((element) => {
    if (names[0].includes(element)) {
      return true;
    }

    return false;
  });
  // console.log(names);
  // console.log(contains);
  let colors = {};

  if (!contains) {
    var randomColors = [];
    while (randomColors.length < names.length) {
      do {
        var color = Math.floor(Math.random() * 10000000000 + 1);
      } while (randomColors.indexOf(color) >= 0);
      randomColors.push('#' + ('000000' + color.toString(16)).slice(-6));
    }
    names.forEach((element, index) => {
      colors[element] = randomColors[index];
    });
  } else {
    colors = {
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
      17: '#ef4c7d',
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
      88: '#BB9139',
      92: '#0E1844',
      110: '#5E1855',
      '-others': '#cececa',
    };
  }

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

  const heatmapColorscale2 = [
    [0, 'rgb(69,8,87)'],
    [0.1, 'rgb(70,36,106)'],

    [0.1, 'rgb(69,44,112)'],
    [0.2, 'rgb(68,55,123)'],

    [0.2, 'rgb(66,64,130)'],
    [0.3, 'rgb(59,92,138)'],

    [0.3, 'rgb(57,97,139)'],
    [0.4, 'rgb(49,111,141)'],

    [0.4, 'rgb(43,119,142)'],
    [0.5, 'rgb(43,139,138)'],

    [0.5, 'rgb(42,151,136'],
    [0.6, 'rgb(35,165,132)'],

    [0.6, 'rgb(39,168,130)'],
    [0.7, 'rgb(84,185,112)'],

    [0.7, 'rgb(101,195,100)'],
    [0.8, 'rgb(120,208,82)'],

    [0.8, 'rgb(145,213,76)'],
    [0.9, 'rgb(177,218,68)'],

    [0.9, 'rgb(204,223,59)'],
    [1.0, 'rgb(249,230,39)'],
  ];

  const heatmapColorscale3 = [
    [0, '#4DBBD5'],
    [0.1, '#3C5488'],
    [0.2, '#E64B35'],

    [0.3, '#F39B7F'],
    [0.4, '#00A087'],

    [0.5, 'EBC6C4'],
    [0.6, 'rgb(35,165,132)'],

    [0.7, 'rgb(101,195,100)'],
    [0.8, 'rgb(120,208,82)'],

    [0.9, 'rgb(204,223,59)'],
    [1.0, 'rgb(249,230,39)'],
  ];

  const xAxis = cosineData.map((e) => e.sample);
  //const longest = xAxis.reduce((a, e) => (a > e.length ? a : e.length), 0);
  var longest = xAxis.sort(function (a, b) {
    return b.length - a.length;
  })[0].length;
  const extraMargin =
    longest > 0 && longest < 10 ? -0.157 : (longest * -0.027) / 2;

  // console.log(xAxis);
  // console.log(longest);
  // console.log(extraMargin);
  function signatureSort(a, b) {
    return a[0].signatureName.localeCompare(b[0].signatureName, 'en', {
      numeric: true,
    });
  }
  let colorBarLoc;
  if (xAxis.length > 250 || arrNames.length > 0) {
    if (longest > 30) {
      colorBarLoc = -0.3;
    } else {
      colorBarLoc = -0.2;
    }
  } else {
    if (longest > 25) {
      colorBarLoc = -0.7;
    } else if (longest > 15) {
      colorBarLoc = -0.4;
    } else {
      colorBarLoc = -0.33;
    }
  }

  console.log(colorBarLoc);
  const groupBySignatureName = groupBy(exposureData, 'signatureName');

  const stackedBarTraces = Object.values(groupBySignatureName)
    .sort(signatureSort)
    .reverse()
    .map((sigArray) => ({
      type: 'bar',
      showlegend: true,
      name: sigArray[0].signatureName,
      x: sigArray.map((e) => e.sample),
      y: sigArray.map((e) => e.exposure),
      marker: {
        color: colors[sigArray[0].signatureName.replace(/^\D*/, '')],
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

  const normalizedBarTraces = Object.values(groupBySignatureName)
    .sort(signatureSort)
    .reverse()
    .map((sigArray) => {
      const groupBySample = groupBy(sigArray, 'sample');
      const totalPerSample = Object.entries(groupBySample).reduce(
        (obj, [sample, data]) => ({
          ...obj,
          [sample]: data.reduce((total, e) => total + e.exposure, 0),
        }),
        {}
      );

      return {
        type: 'bar',
        name: sigArray[0].signatureName,
        x: sigArray.map((e) => e.sample),
        y: sigArray.map((e) => e.exposure / totalPerSample[e.sample]),
        marker: {
          color: colors[sigArray[0].signatureName.replace(/^\D*/, '')],
        },
      };
    });

  // console.log(normalizedBarTraces);

  const sortSignatureName = (sourceArray) => {
    const sortByLocation = (a, b) =>
      a.signatureName.localeCompare(b.signatureName, 'en', { numeric: true });
    return sourceArray.sort(sortByLocation);
  };

  const groupBySample_exposure = groupBy(exposureData, 'sample');
  // console.log(groupBySample_exposure);

  // const sortedCosin = cosineData.sort(
  //   (a, b) => xAxis.indexOf(a.sample) - xAxis.indexOf(b.sample)
  // );
  // console.log(sortedCosin);

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

  // console.log(dataSignature2);

  const groupBySignatureName_exposure2 = groupBy(
    sortSignatureName(dataSignature2),
    //dataSignature2,
    'signatureName'
  );

  const tracesNormalize = Object.entries(groupBySignatureName_exposure2)
    .reverse()
    .map(([key, value]) => ({
      key: key,
      value: value,
      name: key,
      type: 'bar',
      x: value.map((e) => e.sample),
      y: value.map((e, i) => e.exposure / e.total),
      marker: {
        color: contains ? colors[key.replace(/^\D*/, '')] : colors[key],
      },
      showlegend: false,
      exposure: value.map((e, i) => e.exposure),
      hovertemplate:
        '<b>Signature contribution: </b>%{y} <br><b>Sample: </b> %{x}',
    }));
  //console.log(tracesNormalize);
  const tracesHeatMap = [
    {
      z: [cosineData.map((e) => e.similarity)],
      x: cosineData.map((e) => e.sample),
      //hoverongaps: false,
      xaxis: 'x',
      yaxis: 'y2',
      type: 'heatmap',
      colorscale: heatmapColorscale,
      zmin: 0,
      zmax: 1,
      colorbar: {
        orientation: 'h',
        x: 0.5,
        y: colorBarLoc,
        thickness: 20,
        bordercolor: 'black',
        tickmode: 'array',
        tickvals: [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1],
      },
      xgap: 0.5,
      hovertemplate:
        '<b>Sample: </b> %{x}<br> <b>Similarity: </b>%{z} <extra></extra>',
    },
  ];
  console.log(tracesHeatMap);
  let tracesStackedBar;
  if (arrNames.length > 0) {
    tracesStackedBar = Object.entries(groupBySignatureName_exposure2)
      .reverse()
      .map(([key, value]) => ({
        key: key,
        value: value,
        name: key,
        type: 'bar',
        x: value.filter((obj) => obj.exposure !== 0).map((e) => e.sample),
        y: value.filter((obj) => obj.exposure !== 0).map((e, i) => e.exposure),
        marker: {
          color: contains ? colors[key.replace(/^\D*/, '')] : colors[key],
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
        legendgroup: 'b',
        legendgrouptitle: {
          text: '\t Mutational Signatures:',
        },
        showlegend: true,
        hovertemplate:
          '<b>Number of mutation: </b>%{y} <br><b>Sample: </b> %{x}',
      }));
  } else {
    tracesStackedBar = Object.entries(groupBySignatureName_exposure2)
      .reverse()
      .map(([key, value]) => ({
        key: key,
        value: value,
        name: key,
        type: 'bar',
        x: value.filter((obj) => obj.exposure !== 0).map((e) => e.sample),
        y: value.filter((obj) => obj.exposure !== 0).map((e, i) => e.exposure),
        // marker: {
        //   color: colors[key.replace(/^\D*/, '')],
        // },
        marker: {
          color: contains ? colors[key.replace(/^\D*/, '')] : colors[key],
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

        showlegend: true,
        hovertemplate:
          '<b>Number of mutation: </b>%{y} <br><b>Sample: </b> %{x}',
      }));
  }

  const tracesHeatMapVariableNum1 = [
    {
      z: [variableData.map((e) => e.value1)],
      x: variableData.map((e) => e.sample),
      //hoverongaps: false,
      xaxis: 'x',
      yaxis: 'y4',
      type: 'heatmap',
      colorscale: heatmapColorscale2,
      zmin: 0,
      zmax: 100,
      colorbar: {
        orientation: 'h',
        x: 0.5,
        y: colorBarLoc - 0.1,
        thickness: 20,
        bordercolor: 'black',
        tickmode: 'array',
        tickvals: [0, 20, 40, 60, 80, 100],
      },
      xgap: 0.5,
      hovertemplate:
        '<b>Sample: </b> %{x}<br> <b>Value: </b>%{z} <extra></extra>',
    },
  ];

  const tracesHeatMapVariableNum2 = [
    {
      z: [variableData.map((e) => e.value2)],
      x: variableData.map((e) => e.sample),
      //hoverongaps: false,
      xaxis: 'x',
      yaxis: 'y5',
      type: 'heatmap',
      colorscale: heatmapColorscale2,
      zmin: 0,
      zmax: 100,
      showscale: mixMatch ? true : false,
      colorbar: {
        orientation: 'h',
        x: 0.5,
        y: colorBarLoc - 0.1,
        thickness: 20,
        bordercolor: 'black',
        tickmode: 'array',
        tickvals: [0, 20, 40, 60, 80, 100],
      },
      xgap: 0.5,
      hovertemplate:
        '<b>Sample: </b> %{x}<br> <b>Value: </b>%{z} <extra></extra>',
    },
  ];
  // var tracesBarMapVariableStr1 = variableData.map((e) => ({
  //   x: [e.sample],
  //   y: [1],
  //   name: e.value1,
  //   customdata: [
  //     {
  //       name: e.value1,
  //     },
  //   ],
  //   xaxis: 'x',
  //   yaxis: 'y4',
  //   type: 'bar',

  //   marker: { color: charColors[e.value1] },
  //   hovertemplate: '<b>Sample: </b>%{x} <br><b>Value: </b> %{customdata.name}',
  // }));

  const tracesBarMapVariableStr1 = {
    x: variableData.map((e) => e.sample),
    y: variableData.map((e) => 1),
    customdata: variableData.map((e) => ({
      name: e.value1,
    })),
    name: 'Variable String',
    xaxis: 'x',
    yaxis: 'y4',
    type: 'bar',
    test: variableData.map((e) => e.value1),
    marker: { color: variableData.map((e) => charColors[e.value1]) },
    // legendgroup: 'heatmap',
    // legendgrouptitle: {
    //   text: '\t Variable',
    // },
    showlegend: false,
    hovertemplate: '<b>Sample: </b>%{x} <br><b>Value: </b> %{customdata.name}',
  };

  const tracesBarMapVariableStr2 = {
    x: variableData.map((e) => e.sample),
    y: variableData.map((e) => 1),
    customdata: variableData.map((e) => ({
      name: e.value2,
    })),
    xaxis: 'x',
    yaxis: 'y4',
    //type: 'bar',
    test: variableData.map((e) => e.value1),
    marker: { color: variableData.map((e) => charColors[e.value1]) },
    // legendgroup: 'heatmap',
    // legendgrouptitle: {
    //   text: '\t Variable',
    // },
    showlegend: false,
    hovertemplate: '<b>Sample: </b>%{x} <br><b>Value: </b> %{customdata.name}',
  };
  const stringLegend = Object.entries(charColors)
    .sort()
    .reverse()
    .map(([key, val], index) => ({
      key: key,
      val: val,
      x: variableData.map((e) => xAxis[0]),
      y: variableData.map((e) => 0),
      xaxis: 'x',
      yaxis: 'y',
      name: key,
      marker: {
        size: 0,
        color: val,
      },
      type: 'bar',
      hoverinfo: 'skip',
      legendgroup: 'a',
      legendgrouptitle: {
        text: '\t Variable Data String:',
      },
    }));
  console.log(stringLegend);
  console.log(tracesBarMapVariableStr1);
  let traces = [];
  if (variableData.length > 0) {
    if (variableData[0].hasOwnProperty('value2')) {
      if (
        !isNumber(variableData[0].value1) &&
        isNumber(variableData[0].value2)
      ) {
        //1st is string , 2nd is numeric
        traces = [
          ...tracesHeatMap,
          ...tracesStackedBar,
          ...tracesNormalize,
          tracesBarMapVariableStr1,
          ...stringLegend,
          ...tracesHeatMapVariableNum2,
        ];
      } else if (
        isNumber(variableData[0].value1) &&
        isNumber(variableData[0].value2)
      ) {
        //both value are number
        traces = [
          ...tracesHeatMap,
          ...tracesStackedBar,
          ...tracesNormalize,
          ...tracesHeatMapVariableNum1,
          ...tracesHeatMapVariableNum2,
        ];
      } else {
        //both value are number
        traces = [
          ...tracesHeatMap,
          ...tracesStackedBar,
          ...tracesNormalize,
          tracesBarMapVariableStr1,
          ...stringLegend,
          tracesBarMapVariableStr2,
        ];
      }
    } else {
      if (!isNumber(variableData[0].value1)) {
        traces = [
          ...tracesHeatMap,
          ...tracesStackedBar,
          ...tracesNormalize,
          tracesBarMapVariableStr1,
          ...stringLegend,
        ];
      } else {
        traces = [
          ...tracesHeatMap,
          ...tracesStackedBar,
          ...tracesNormalize,
          ...tracesHeatMapVariableNum1,
        ];
      }
    }
  } else {
    traces = [...tracesHeatMap, ...tracesStackedBar, ...tracesNormalize];
  }

  const text = {
    x: 0,
    y: xAxis.length > 250 ? -0.03 : extraMargin,

    xanchor: 'left',
    yanchor: 'bottom',
    xref: 'paper',
    yref: 'paper',
    text: '\t Mutational Signatures:',
    showarrow: false,
    font: {
      family: 'Arial',
      size: 17,
      color: 'rgb(37,37,37)',
    },
  };

  const annotationTitle = {
    x: 0,
    y: colorBarLoc + 0.06,

    xanchor: 'left',
    yanchor: 'bottom',
    xref: 'paper',
    yref: 'paper',
    text: '\t Cosine Similarity:',
    showarrow: false,
    font: {
      family: 'Arial',
      size: 17,
      color: 'rgb(37,37,37)',
    },
  };

  const annotationLegendTitle = {
    x: 0,
    y: colorBarLoc - 0.043,

    xanchor: 'left',
    yanchor: 'bottom',
    xref: 'paper',
    yref: 'paper',
    text: '\t Variable Data Number:',
    showarrow: false,
    font: {
      family: 'Arial',
      size: 17,
      color: 'rgb(37,37,37)',
    },
  };

  let annotations =
    variableData.length > 0
      ? arrNames.length > 0
        ? mixMatch
          ? [annotationTitle, annotationLegendTitle]
          : [annotationTitle]
        : [text, annotationTitle, annotationLegendTitle]
      : [text, annotationTitle];

  const layout = {
    title: {
      text: '<b>Landscape of Mutational Signature Activity</b>',
      font: {
        family: 'Arial',
        size: 18,
      },
    },
    autosize: true,
    height: 1200,
    barmode: 'stack',
    hovermode: 'closest',
    // legend: {
    //   orientation: 'h',
    //   //title: { text: 'Signatures Name<br>' },
    //   traceorder: 'reserved',
    //   x: 0,
    //   y: xAxis.length > 250 ? -0.03 : extraMargin,
    // },
    legend: {
      orientation: arrNames.length > 0 ? 'v' : 'h',
      x: arrNames.length > 0 ? 1 : 0,
      y: arrNames.length > 0 ? 1 : xAxis.length > 250 ? -0.03 : extraMargin,
    },

    xaxis: {
      tickmode: 'array',
      tickvals: xAxis.map((e) => e),
      // tickvals: xAxis.map((_, i) => i),
      ticktext: xAxis.map((e) => e),
      type: 'category',
      tickangle: -90,
      ticks: '',
      showticklabels: xAxis.length > 230 ? false : true,
      showgrid: false,
      zeroline: false,
      showline: false,
      zerolinewidth: 0,
    },
    yaxis: { title: 'Signature contribution', domain: [0, 0.49] },
    yaxis2: {
      title: '',
      domain: [0.475, 0.493],
      showticklabels: false,
      ticks: '',
    },
    yaxis3: { title: 'Number of mutation', domain: [0.5, 0.95] },
    yaxis4: {
      title: '',
      domain: [0.96, 0.978],
      showticklabels: false,
      ticks: '',
    },
    yaxis5: {
      title: '',
      domain: [0.978, 0.996],
      showticklabels: false,
      ticks: '',
    },

    bargap: 0.04,
    annotations: annotations,
  };

  var config = {
    //responsive: true,
  };

  return { traces, layout, config };
}
