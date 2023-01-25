import { groupBy } from 'lodash';
import { First } from 'react-bootstrap/esm/PageItem';
import { get } from 'react-hook-form';
import CosineSimilarityPlot from '../../../pages/catalog/referenceSignature/cosineSimilarity/form-plot';
import { getRandomColor, isNumber, mapOrder } from '../../utils/utils';
import { colorPallet, colorPallet0 } from '../../utils/colors';

export default function MsLandscape(cosineData, exposureData, variableData) {
  console.log(cosineData);
  console.log(exposureData);
  console.log(variableData);

  let charColors = {};
  let arrNames = [];
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

      stringData.forEach((element, index) => {
        charColors[element] = bg_colors[index];
      });
      arrNames = stringData;
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

  const defaultNames = ['SBS', 'DBS', 'ID'];
  const names = exposureData.map((group) => group.signatureName);
  const contains = defaultNames.some((element) => {
    if (names[0].includes(element)) {
      return true;
    }

    return false;
  });
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
    colors = colorPallet;
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

  const newVariableData = [...variableData];
  const variableDataSort = newVariableData.sort(mapOrder(xAxis, 'sample'));

  // console.log(longest);
  // console.log(extraMargin);
  function signatureSort(a, b) {
    return a[0].signatureName.localeCompare(b[0].signatureName, 'en', {
      numeric: true,
    });
  }
  let colorBarLoc;
  // if (xAxis.length > 250 || arrNames.length > 0) {
  //   if (longest > 30) {
  //     colorBarLoc = -0.3;
  //   } else {
  //     colorBarLoc = -0.2;
  //   }
  // } else {
  //   if (longest > 25) {
  //     colorBarLoc = -0.7;
  //   } else if (longest > 15) {
  //     colorBarLoc = -0.4;
  //   } else {
  //     colorBarLoc = -0.33;
  //   }
  // }
  if (longest > 200) {
    colorBarLoc = -0.2;
  } else if (longest > 30) {
    colorBarLoc = -0.3;
  } else {
    colorBarLoc = -0.2;
  }

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
        legendrank: 1001,
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
        legendgrouptitle: {
          text: '\t Mutational Signatures:',
        },
        showlegend: true,
        hovertemplate:
          '<b>Number of mutation: </b>%{y} <br><b>Sample: </b> %{x}',
      }));
  }

  const tracesHeatMapVariableNum1 = [
    {
      z: [variableDataSort.map((e) => e.value1)],
      x: variableDataSort.map((e) => e.sample),
      //hoverongaps: false,
      xaxis: variableDataSort.length < xAxis.length ? 'x2' : 'x',
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
      xAxis: xAxis.length,
      variableData: variableDataSort.length,
      xgap:
        variableDataSort.length < xAxis.length
          ? Math.ceil(xAxis.length / variableDataSort.length) * 2.5
          : 0.5,
      //xgap: 2,
      hovertemplate:
        '<b>Sample: </b> %{x}<br> <b>Value: </b>%{z} <extra></extra>',
    },
  ];
  console.log(tracesHeatMapVariableNum1);
  const tracesHeatMapVariableNum2 = [
    {
      z: [variableDataSort.map((e) => e.value2)],
      x: variableDataSort.map((e) => e.sample),
      //hoverongaps: false,
      xaxis: variableDataSort.length < xAxis.length ? 'x2' : 'x',
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
      xgap:
        variableDataSort.length < xAxis.length
          ? Math.ceil(xAxis.length / variableDataSort.length) * 2.5
          : 0.5,
      //xgap: 2,
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
    x: variableDataSort.map((e) => e.sample),
    y: variableDataSort.map((e) => 1),
    customdata: variableDataSort.map((e) => ({
      name: e.value1,
      xVal: e.sample,
    })),

    xaxis: variableDataSort.length < xAxis.length ? 'x2' : 'x',
    yaxis: 'y4',
    type: 'bar',
    test: variableDataSort.map((e) => e.value1),
    marker: {
      color: variableDataSort.map((e) => charColors[e.value1]),
      //line: { width: 1 },
    },
    showlegend: false,
    hovertemplate:
      '<b>Sample: </b>%{x} <br><b>Value: </b> %{customdata.name}<extra></extra>',
  };

  const tracesBarMapVariableStr2 = {
    x: variableDataSort.map((e) => e.sample),
    y: variableDataSort.map((e) => 1),
    customdata: variableDataSort.map((e) => ({
      name: e.value2,
    })),
    xaxis: variableDataSort.length < xAxis.length ? 'x2' : 'x',
    yaxis: 'y4',
    type: 'bar',
    test: variableDataSort.map((e) => e.value1),
    marker: {
      color: variableDataSort.map((e) => charColors[e.value1]),
      //line: { width: 1 },
    },
    showlegend: false,
    hovertemplate:
      '<b>Sample: </b>%{x} <br><b>Value: </b> %{customdata.name}<extra></extra>',
  };
  const stringLegend = Object.entries(charColors)
    .sort()
    .reverse()
    .map(([key, val], index) => ({
      key: key,
      val: val,
      x: variableDataSort.map((e) => xAxis[0]),
      y: variableDataSort.map((e) => 0),
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

  let traces = [];
  if (variableDataSort.length > 0) {
    if (variableDataSort[0].hasOwnProperty('value2')) {
      if (
        !isNumber(variableDataSort[0].value1) &&
        isNumber(variableDataSort[0].value2)
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

  // let annotations =
  //   variableData.length > 0
  //     ? arrNames.length > 0
  //       ? mixMatch
  //         ? [annotationTitle, annotationLegendTitle]
  //         : [annotationTitle]
  //       : [text, annotationTitle, annotationLegendTitle]
  //     : [text, annotationTitle];

  let annotations =
    variableData.length > 0
      ? arrNames.length > 0
        ? mixMatch
          ? [annotationTitle, annotationLegendTitle]
          : [annotationTitle]
        : [annotationTitle, annotationLegendTitle]
      : [annotationTitle];
  const layout = {
    title: {
      text: '<b>Landscape of Mutational Signature Activity</b>',
      font: {
        family: 'Arial',
        size: 18,
      },
    },
    bargap: 0.1,
    autosize: true,
    height: 1200,
    barmode: 'stack',
    hovermode: 'closest',
    legend: {
      x: 1,
      y: 0.95,
    },
    // legend: {
    //   orientation: arrNames.length > 0 ? 'v' : 'h',
    //   x: arrNames.length > 0 ? 1 : 0,
    //   y: arrNames.length > 0 ? 1 : xAxis.length > 250 ? -0.03 : extraMargin,
    // },

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
    xaxis2: {
      showticklabels: false,
      zerolinecolor: 'rgba(0,0,0,0)',
    },
    yaxis: {
      title: 'Signature contribution',
      domain: [0, 0.49],
      zeroline: false,
    },
    yaxis2: {
      title: '',
      domain: [0.475, 0.493],
      showticklabels: false,
      ticks: '',
      showgrid: true,
      zeroline: false,
    },
    yaxis3: {
      title: 'Number of mutation',
      domain: [0.5, 0.95],
      zeroline: false,
    },
    yaxis4: {
      title: '',
      domain: [0.96, 0.978],
      showticklabels: false,
      ticks: '',
      zeroline: false,
    },
    yaxis5: {
      title: '',
      domain: [0.978, 0.996],
      showticklabels: false,
      ticks: '',
      zeroline: false,
    },

    annotations: annotations,
  };

  var config = {
    //responsive: true,
  };

  return { traces, layout, config };
}
