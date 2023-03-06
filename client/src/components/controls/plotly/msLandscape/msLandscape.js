import { groupBy } from 'lodash';
import { isNumber, mapOrder } from '../../utils/utils';
import { colorPallet, colorPallet0 } from '../../utils/colors';

export default function MsLandscape(
  cosineData,
  exposureData,
  variableData,
  dendrogram = {}
) {
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

  const samples = cosineData.map((e) => e.sample);

  var longest = samples.sort(function (a, b) {
    return b.length - a.length;
  })[0].length;

  const newVariableData = [...variableData];
  const variableDataSort = newVariableData.sort(mapOrder(samples, 'sample'));

  let colorBarLoc;
  if (longest > 200) {
    colorBarLoc = -0.2;
  } else if (longest > 30) {
    colorBarLoc = -0.3;
  } else {
    colorBarLoc = -0.2;
  }

  const sortSignatureName = (sourceArray) => {
    const sortByLocation = (a, b) =>
      a.signatureName.localeCompare(b.signatureName, 'en', { numeric: true });
    return sourceArray.sort(sortByLocation);
  };

  const groupBySample_exposure = groupBy(exposureData, 'sample');

  const dataSignature2 = Object.entries(groupBySample_exposure)
    .map(([sample, data]) => {
      const total = data.reduce((a, e) => a + parseInt(e.exposure), 0);

      return data.map((e) => ({ ...e, total }));
    })
    .flat();

  const groupBySignatureName_exposure2 = groupBy(
    sortSignatureName(dataSignature2),
    'signatureName'
  );

  const signatureContributionTraces = Object.entries(
    groupBySignatureName_exposure2
  )
    .reverse()
    .map(([signature, data]) => ({
      name: signature,
      type: 'bar',
      x: data.map((e) => samples.findIndex((v) => v === e.sample)),
      y: data.map((e, i) => e.exposure / e.total),
      marker: {
        color: contains
          ? colors[signature.replace(/^\D*/, '')]
          : colors[signature],
      },
      showlegend: false,
      exposure: data.map((e, i) => e.exposure),
      hovertemplate:
        '<b>Signature Contribution: </b>%{y}<br><b>Sample: </b>%{x}',
    }));

  const cosineSimHeatMap = [
    {
      z: [cosineData.map((e) => e.similarity)],
      x: cosineData.map((e) => samples.findIndex((v) => v === e.sample)),
      customdata: [
        cosineData.map((e) => ({
          sample: e.sample,
        })),
      ],
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
        '<b>Sample: </b>%{customdata.sample}<br><b>Cosine Similarity: </b>%{z}<extra></extra>',
    },
  ];
  const mutationCountTraces = Object.entries(groupBySignatureName_exposure2)
    .reverse()
    .map(([key, value]) => ({
      name: key,
      type: 'bar',
      x: value
        .filter((obj) => obj.exposure !== 0)
        .map((e) => samples.findIndex((v) => v === e.sample)),
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
      ...(arrNames.length > 0 && {
        legendrank: 1001,
        legendgroup: 'b',
      }),
      legendgrouptitle: {
        text: '\t Mutational Signatures:',
      },
      showlegend: true,
      hovertemplate: '<b>Number of mutation: </b>%{y} <br><b>Sample: </b>%{x}',
    }));

  const tracesHeatMapVariableNum1 = [
    {
      z: [variableDataSort.map((e) => e.value1)],
      x: variableDataSort.map((e, i) => i),
      customdata: [
        variableDataSort.map((e) => ({
          sample: e.sample,
          value: e.value1,
        })),
      ],
      xaxis: variableDataSort.length === samples.length ? 'x' : 'x2',
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
      variableData: variableDataSort.length,
      xgap:
        variableDataSort.length < samples.length
          ? Math.ceil(samples.length / variableDataSort.length) * 2.5
          : 0.5,
      hovertemplate:
        '<b>Sample: </b>%{customdata.sample}<br><b>Value: </b>%{z}<extra></extra>',
    },
  ];

  const tracesHeatMapVariableNum2 = [
    {
      z: [variableDataSort.map((e) => e.value2)],
      x: variableDataSort.map((e, i) => i),
      customdata: [
        variableDataSort.map((e) => ({
          sample: e.sample,
          value: e.value1,
        })),
      ],
      xaxis: variableDataSort.length === samples.length ? 'x' : 'x2',
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
        variableDataSort.length < samples.length
          ? Math.ceil(samples.length / variableDataSort.length) * 2.5
          : 0.5,
      hovertemplate:
        '<b>Sample: </b>%{customdata.sample}<br><b>Value: </b>%{z}<extra></extra>',
    },
  ];

  const tracesBarMapVariableStr1 = {
    x: variableDataSort.map((e, i) => i),
    y: variableDataSort.map((e) => 1),
    customdata: variableDataSort.map((e) => ({
      sample: e.sample,
      value: e.value1,
    })),
    xaxis: variableDataSort.length === samples.length ? 'x' : 'x2',
    yaxis: 'y4',
    type: 'bar',
    marker: {
      color: variableDataSort.map((e) => charColors[e.value1]),
    },
    showlegend: false,
    hovertemplate:
      '<b>Sample: </b>%{customdata.sample}<br><b>Value: </b>%{customdata.value}<extra></extra>',
  };

  const tracesBarMapVariableStr2 = {
    x: variableDataSort.map((e, i) => i),
    y: variableDataSort.map((e) => 1),
    customdata: variableDataSort.map((e) => ({
      name: e.value2,
    })),
    xaxis: variableDataSort.length === samples.length ? 'x' : 'x2',
    yaxis: 'y4',
    type: 'bar',
    marker: {
      color: variableDataSort.map((e) => charColors[e.value1]),
    },
    showlegend: false,
    hovertemplate:
      '<b>Sample: </b>%{x}<br><b>Value: </b>%{customdata.name}<extra></extra>',
  };
  const stringLegend = Object.entries(charColors)
    .sort()
    .reverse()
    .map(([key, val], index) => ({
      x: variableDataSort.map((e) => samples[0]),
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
        text:
          variableData.length > 0
            ? '\t' + variableData[0].titles[0]
            : '\t Variable Data String:',
      },
    }));

  // modify dendrogram trace
  const dendrogramTrace = Object.keys(dendrogram).length
    ? {
        ...dendrogram.data[1],
        xaxis: 'x3',
        yaxis: 'y6',
        text: dendrogram.data[1].x.map((e) =>
          samples[e] ? `<b>Sample: </b>${samples[e - 1]}` : null
        ),
      }
    : {};

  let traces = [];
  if (variableDataSort.length > 0) {
    if (variableDataSort[0].hasOwnProperty('value2')) {
      if (
        !isNumber(variableDataSort[0].value1) &&
        isNumber(variableDataSort[0].value2)
      ) {
        //1st is string , 2nd is numeric
        traces = [
          ...cosineSimHeatMap,
          ...mutationCountTraces,
          ...signatureContributionTraces,
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
          ...cosineSimHeatMap,
          ...mutationCountTraces,
          ...signatureContributionTraces,
          ...tracesHeatMapVariableNum1,
          ...tracesHeatMapVariableNum2,
        ];
      } else {
        //both value are number
        traces = [
          ...cosineSimHeatMap,
          ...mutationCountTraces,
          ...signatureContributionTraces,
          tracesBarMapVariableStr1,
          ...stringLegend,
          tracesBarMapVariableStr2,
        ];
      }
    } else {
      if (!isNumber(variableData[0].value1)) {
        traces = [
          ...cosineSimHeatMap,
          ...mutationCountTraces,
          ...signatureContributionTraces,
          tracesBarMapVariableStr1,
          ...stringLegend,
        ];
      } else {
        traces = [
          ...cosineSimHeatMap,
          ...mutationCountTraces,
          ...signatureContributionTraces,
          ...tracesHeatMapVariableNum1,
        ];
      }
    }
  } else {
    traces = [
      ...cosineSimHeatMap,
      ...mutationCountTraces,
      ...signatureContributionTraces,
    ];
  }

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
  let variableDataNum = '';
  if (variableData.length > 0) {
    if (variableData[0].titles.length === 2) {
      if (arrNames.length > 0) {
        if (mixMatch) {
          variableDataNum = '\t' + variableData[0].titles[1];
        } else {
          variableDataNum = '\t' + variableData[0].titles[0];
        }
      } else {
        variableDataNum =
          '\t' + variableData[0].titles[0] + '/ ' + variableData[0].titles[1];
      }
    } else {
      variableDataNum = '\t' + variableData[0].titles[0];
    }
  }

  const annotationLegendTitle = {
    x: 0,
    y: colorBarLoc - 0.043,

    xanchor: 'left',
    yanchor: 'bottom',
    xref: 'paper',
    yref: 'paper',
    text: variableDataNum,
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
    height: 1500,
    barmode: 'stack',
    hovermode: 'closest',
    legend: {
      x: 1,
      y: 0.95,
    },

    xaxis: {
      tickmode: 'array',
      tickvals: samples.map((e, i) => i),
      ticktext: samples,
      type: 'category',
      tickangle: -90,
      ticks: '',
      showticklabels: samples.length > 230 ? false : true,
      showgrid: false,
      zeroline: false,
      showline: false,
      zerolinewidth: 0,
    },
    xaxis2: {
      showticklabels: false,
      zerolinecolor: 'rgba(0,0,0,0)',
    },
    xaxis3: {
      ...(Object.keys(dendrogram).length && dendrogram.layout.xaxis),
      showticklabels: false,
      range: [0.5, samples.length + 0.5],
    },
    yaxis: {
      title: 'Signature Contribution',
      domain: [0, 0.35],
      zeroline: false,
    },
    yaxis2: {
      domain: [0.355, 0.38],
      showticklabels: false,
      ticks: '',
      showgrid: true,
      zeroline: false,
    },
    yaxis3: {
      title: 'Number of mutation',
      domain: [0.395, 0.74],
      zeroline: false,
    },
    yaxis4: {
      domain: [0.75, 0.77],
      showticklabels: false,
      ticks: '',
      zeroline: false,
    },
    yaxis5: {
      domain: [0.77, 0.79],
      showticklabels: false,
      ticks: '',
      zeroline: false,
    },
    yaxis6: {
      ...(Object.keys(dendrogram).length && dendrogram.layout.yaxis),
      domain: [0.8, 1],
      zeroline: false,
    },

    annotations: annotations,
  };

  const config = { responsive: true };

  return {
    traces: [...traces, dendrogramTrace],
    layout,
    config,
  };
}
