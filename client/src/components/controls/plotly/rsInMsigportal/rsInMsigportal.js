import { groupBy } from 'lodash';
export default function RsInMsigportal(rawData) {
  console.log(rawData);
  function mapOrder(array, order, key) {
    array.sort(function (a, b) {
      var A = a[key],
        B = b[key];

      if (order.indexOf(A) > order.indexOf(B)) {
        return 1;
      } else {
        return -1;
      }
    });

    return array;
  }
  const profile_order = ['SBS', 'DBS', 'ID', 'CN', 'RS'];
  const groupBySpecies = groupBy(rawData, (item) => `${item.species}`);
  const groupByignatureSetName = groupBy(
    rawData,
    (item) => `${item.signatureSetName}`
  );

  const signatureSetName = Object.keys(groupByignatureSetName).map((e) => e);

  var randomColors = [];
  while (randomColors.length < signatureSetName.length) {
    do {
      var color = Math.floor(Math.random() * 10000000000 + 1);
    } while (randomColors.indexOf(color) >= 0);
    randomColors.push('#' + ('000000' + color.toString(16)).slice(-6));
  }
  let colors = {};
  Object.keys(groupByignatureSetName).map((e, i) => {
    colors[e] = randomColors[i];
  });
  let dataHuman = [];
  let dataMm9 = [];
  let dataRn6 = [];
  Object.entries(groupBySpecies).map(([key, val], index) => {
    if (key.includes('(GRCh37/38)')) {
      dataHuman.push(val);
    } else if (key.includes('(mm9/10)')) {
      dataMm9.push(val);
    } else {
      dataRn6.push(val);
    }
  });
  dataHuman = dataHuman.flat();
  const dataHumanSort = mapOrder(dataHuman, profile_order, 'profile');
  dataMm9 = dataMm9.flat();
  dataRn6 = dataRn6.flat();

  console.log(dataHuman);
  console.log(dataHumanSort);

  const groupedHuman = groupBy(
    dataHuman,
    (item) => `${item.profile}${item.matrix}`
  );
  console.log(groupedHuman);
  const groupedMm9 = groupBy(
    dataMm9,
    (item) => `${item.profile}${item.matrix}`
  );
  const groupedRn6 = groupBy(
    dataRn6,
    (item) => `${item.profile}${item.matrix}`
  );
  // const dataArray = [];
  // Object.values(groupedHuman).map((e) => dataArray.push(e));

  for (var key in groupedHuman)
    groupedHuman[key].sort(function (a, b) {
      var x = a.profile,
        y = b.profile;
      return x === y ? 0 : x < y ? -1 : 1;
    });
  console.log(groupedHuman);

  const tracePies0 = Object.entries(groupedHuman).map(
    ([key, element], index, array) => ({
      type: 'pie',
      marker: {
        color: element.map((e) => colors[e.signatureSetName]),
        line: {
          color: 'black',
          width: 1,
        },
      },
      textposition: 'inside',
      labels: element.map((e) => e.signatureSetName),
      values: element.map((e) => parseInt(e.count)),
      texttemplate: '%{value}',
      direction: 'clockwise',
      name: key,
      domain: {
        x: [
          index < 4
            ? Math.round(index * (1 / 5) * 10) / 10
            : Math.round((index - 4) * (1 / 5) * 10) / 10,
          index < 4
            ? Math.round((index * (1 / 5) + 0.2) * 10) / 10
            : Math.round(((index - 4) * (1 / 5) + 0.2) * 10) / 10,
        ],
        y: [index < 4 ? 0.55 : 0.78, index < 4 ? 0.745 : 0.975],
      },
      hovertemplate:
        '<b>%{label}</b> <br>%{percent} </br> %{value}  <extra></extra>',
    })
  );
  const tracePies1 = Object.entries(groupedMm9).map(
    ([key, element], index, array) => ({
      type: 'pie',
      e: element,
      marker: {
        color: element.map((e) => colors[e.signatureSetName]),
        line: {
          color: 'black',
          width: 1,
        },
      },
      textposition: 'inside',
      labels: element.map((e) => e.signatureSetName),
      values: element.map((e) => parseInt(e.count)),
      texttemplate: '%{value}',
      direction: 'clockwise',
      name: key,
      domain: {
        x: [
          index < 4
            ? Math.round(index * (1 / 5) * 10) / 10
            : Math.round((index - 4) * (1 / 5) * 10) / 10,
          index < 4
            ? Math.round((index * (1 / 5) + 0.2) * 10) / 10
            : Math.round(((index - 4) * (1 / 5) + 0.2) * 10) / 10,
        ],
        y: [0.265, 0.46],
      },
      hovertemplate:
        '<b>%{label}</b> <br>%{percent} </br> %{value}  <extra></extra>',
    })
  );

  const tracePies2 = Object.entries(groupedRn6).map(
    ([key, element], index, array) => ({
      type: 'pie',
      e: element,
      marker: {
        color: element.map((e) => colors[e.signatureSetName]),
        line: {
          color: 'black',
          width: 1,
        },
      },
      textposition: 'inside',
      labels: element.map((e) => e.signatureSetName),
      values: element.map((e) => parseInt(e.count)),
      texttemplate: '%{value}',
      direction: 'clockwise',
      name: key,
      domain: {
        x: [
          index < 4
            ? Math.round(index * (1 / 5) * 10) / 10
            : Math.round((index - 4) * (1 / 5) * 10) / 10,
          index < 4
            ? Math.round((index * (1 / 5) + 0.2) * 10) / 10
            : Math.round(((index - 4) * (1 / 5) + 0.2) * 10) / 10,
        ],
        y: [0, 0.195],
      },
      hovertemplate:
        '<b>%{label}</b> <br>%{percent} </br> %{value}  <extra></extra>',
    })
  );
  function indexPos(index) {
    let indexPosition;
    if (index === 0 || index === 4) {
      indexPosition = 0.072;
    } else if (index === 1 || index === 5) {
      indexPosition = 0.27;
    } else if (index === 2 || index === 6) {
      indexPosition = 0.495;
    } else if (index === 3 || index === 7) {
      indexPosition = 0.72;
    } else {
      indexPosition = 0.92;
    }
    return indexPosition;
  }

  const pieTitles0 = Object.entries(groupedHuman).map(
    ([key, element], index, array) => ({
      xref: 'paper',
      yref: 'paper',
      xanchor: 'bottom',
      yanchor: 'bottom',
      showarrow: false,
      text: key.padStart(7, ' '),
      align: 'center',
      //x: index < 4 ? index * (1 / 5) + 0.0975 : (index - 4) * (1 / 5) + 0.0975,
      x: indexPos(index),
      y: index < 4 ? 0.75 : 0.98,
    })
  );
  console.log(pieTitles0);
  const pieTitles1 = Object.entries(groupedMm9).map(
    ([key, element], index, array) => ({
      xref: 'paper',
      yref: 'paper',
      xanchor: 'bottom',
      yanchor: 'bottom',
      showarrow: false,
      text: key.padStart(7, ' '),
      align: 'right',
      // x: index < 4 ? index * (1 / 5) + 0.085 : (index - 4) * (1 / 5) + 0.085,
      x: indexPos(index),
      y: 0.46,
    })
  );

  const pieTitles2 = Object.entries(groupedRn6).map(
    ([key, element], index, array) => ({
      xref: 'paper',
      yref: 'paper',
      xanchor: 'bottom',
      yanchor: 'bottom',
      showarrow: false,
      text: key.padStart(7, ' '),
      align: 'center',
      //x: index < 4 ? index * (1 / 5) + 0.08 : (index - 4) * (1 / 5) + 0.08,
      x: indexPos(index),
      y: 0.2,
    })
  );
  const annotationTitle0 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    showarrow: false,
    text: Object.keys(groupBySpecies)[0],
    font: {
      size: 16,
    },
    x: 0.5,
    y: 1.02,
  };

  const annotationTitle1 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    showarrow: false,
    text: Object.keys(groupBySpecies)[1],
    font: {
      size: 16,
    },
    x: 0.25,
    y: 0.5,
  };

  const annotationTitle2 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    showarrow: false,
    text: Object.keys(groupBySpecies)[2],
    font: {
      size: 16,
    },
    x: 0.15,
    y: 0.23,
  };

  const shapes = [
    {
      type: 'line',
      xref: 'paper',
      yref: 'paper',
      x0: 0,
      y0: 1.03,
      x1: 0.35,
      y1: 1.03,
      line: {
        color: 'gray',
        width: 3,
      },
    },
    {
      type: 'line',
      xref: 'paper',
      yref: 'paper',
      x0: 0.64,
      y0: 1.03,
      x1: 1,
      y1: 1.03,
      line: {
        color: 'gray',
        width: 3,
      },
    },
    {
      type: 'line',
      xref: 'paper',
      yref: 'paper',
      x0: 0,
      y0: 0.51,
      x1: 0.23,
      y1: 0.51,
      line: {
        color: 'gray',
        width: 3,
      },
    },
    {
      type: 'line',
      xref: 'paper',
      yref: 'paper',
      x0: 0.5,
      y0: 0.51,
      x1: 0.7,
      y1: 0.51,
      line: {
        color: 'gray',
        width: 3,
      },
    },
    {
      type: 'line',
      xref: 'paper',
      yref: 'paper',
      x0: 0,
      y0: 0.24,
      x1: 0.14,
      y1: 0.24,
      line: {
        color: 'gray',
        width: 3,
      },
    },
    {
      type: 'line',
      xref: 'paper',
      yref: 'paper',
      x0: 0.38,
      y0: 0.24,
      x1: 0.7,
      y1: 0.24,
      line: {
        color: 'gray',
        width: 3,
      },
    },
  ];

  const traces = [...tracePies0, ...tracePies1, ...tracePies2];
  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: 975,
    width: 1080,
    autosize: true,
    legend: {
      title: { text: '\t Signature Set Name' },
      x: 1,
      xanchor: 'right',
      y: 0,
    },
    annotations: [
      ...pieTitles0,
      ...pieTitles1,
      ...pieTitles2,
      annotationTitle0,
      annotationTitle1,
      annotationTitle2,
    ],
    shapes: [...shapes],
  };
  const config = {
    //responsive: true,
  };
  return { traces: traces, layout: layout, config };
}
