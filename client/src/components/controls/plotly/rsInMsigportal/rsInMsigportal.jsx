import { groupBy } from 'lodash';
export default function RsInMsigportal(rawData) {
  const profile_order_0 = [
    'SBS96',
    'SBS192',
    'SBS288',
    'SBS1536',
    'DBS78',
    'ID29',
    'ID83',
    'CN48',
    'RS32',
    'RNA192',
  ];
  const profile_order_1 = ['SBS96', 'DBS78', 'ID83'];
  const groupBySpecies = groupBy(rawData, (item) => `${item.species}`);
  console.log("groupBySpecies ", groupBySpecies)
  console.log("rawData in RsInMsigportal: ", rawData);
  const groupByignatureSetName = groupBy(
    rawData,
    (item) => `${item.signatureSetName}`
  );

  console.log("groupByignatureSetName ", groupByignatureSetName);

  const signatureSetName = Object.keys(groupByignatureSetName).map((e) => e);
  console.log("signatureSetName ", signatureSetName);

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
  let dataGallus = [];
  let dataCaenorhabditis = [];
  Object.entries(groupBySpecies).map(([key, val], index) => {
    if (key.includes('(GRCh37/38)')) {
      dataHuman.push(val);
    } else if (key.includes('(mm9/10)')) {
      dataMm9.push(val);
    } else if (key.includes('(rn6)')) {
      dataRn6.push(val);
    } else if (key.includes('GRCg7b')){  
       dataGallus.push(val);
    } else {
      dataCaenorhabditis.push(val);
    }
  });

  console.log("groupBySpecies ", groupBySpecies);
  console.log("dataHuman ", dataHuman);
  console.log("dataMm9 ", dataMm9);
  console.log("dataRn6 ", dataRn6);
  console.log("dataGallus ", dataGallus);

  dataHuman = dataHuman.flat();
  dataMm9 = dataMm9.flat();
  dataRn6 = dataRn6.flat();
  dataGallus = dataGallus.flat();
  dataCaenorhabditis = dataCaenorhabditis.flat();

  const groupedHuman = groupBy(
    dataHuman,
    (item) => `${item.profile}${item.matrix}`
  );
  const groupedMm9 = groupBy(
    dataMm9,
    (item) => `${item.profile}${item.matrix}`
  );
  const groupedRn6 = groupBy(
    dataRn6,
    (item) => `${item.profile}${item.matrix}`
  );
  const groupedGallus = groupBy(
    dataGallus,
    (item) => `${item.profile}${item.matrix}`
  );
  const groupCaenorhabditis = groupBy(
    dataCaenorhabditis,
    (item) => `${item.profile}${item.matrix}`
  );

  const groupedHuman_sorted = profile_order_0.map((key) => groupedHuman[key]);
  console.log("groupedHuman_sorted ", groupedHuman_sorted);
  const groupedMm9_sorted = profile_order_1.map((key) => groupedMm9[key]);
  const tracePies0 = Object.entries(groupedHuman_sorted).map(
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
      index: index,
      name: array[index][1][0].profile + array[index][1][0].matrix,
      domain: {
        x: [
          index < 4 ? 0.02 : index < 8 ? 0.158 : 0.296,
          index < 4 ? 0.116 : index < 8 ? 0.254 : 0.392,
        ],
        y: [
          index < 4 
            ? 1 - ((index + 1) * 0.25) + 0.02
            : index < 8
            ? 1 - ((index - 3) * 0.25) + 0.02
            : 1 - ((index - 7) * 0.25) + 0.02,
          index < 4 
            ? 1 - (index * 0.25) - 0.02
            : index < 8
            ? 1 - ((index - 4) * 0.25) - 0.02
            : 1 - ((index - 8) * 0.25) - 0.02,
        ],
      },
      hovertemplate:
        '<b>%{label}</b> <br>%{percent} </br> %{value}  <extra></extra>',
    })
  );
  console.log("tracePies0 ", tracePies0);
  const tracePies1 = Object.entries(groupedMm9_sorted).map(
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
        x: [0.434, 0.53],
        y: [
          1 - ((index + 1) * 0.25) + 0.02,
          1 - (index * 0.25) - 0.02,
        ],
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
        x: [0.572, 0.668],
        y: [
          1 - ((index + 1) * 0.25) + 0.02,
          1 - (index * 0.25) - 0.02,
        ],
      },
      hovertemplate:
        '<b>%{label}</b> <br>%{percent} </br> %{value}  <extra></extra>',
    })
  );
  const tracePies3 = Object.entries(groupedGallus).map(
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
        x: [0.710, 0.806],
        y: [
          1 - ((index + 1) * 0.25) + 0.02,
          1 - (index * 0.25) - 0.02,
        ],
      },
      hovertemplate:
        '<b>%{label}</b> <br>%{percent} </br> %{value}  <extra></extra>',
    })
  );

  const tracePies4 = Object.entries(groupCaenorhabditis).map(
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
        x: [0.848, 0.944],
        y: [
          1 - ((index + 1) * 0.25) + 0.02,
          1 - (index * 0.25) - 0.02,
        ],
      },
      hovertemplate:
        '<b>%{label}</b> <br>%{percent} </br> %{value}  <extra></extra>',
    })
  );

  function indexPos(index) {
    let indexPosition;
    if (index === 0) {
      indexPosition = 1/12; // Center of first column (0 to 1/6)
    } else if (index === 1) {
      indexPosition = 3/12; // Center of second column (1/6 to 2/6)
    } else if (index === 2) {
      indexPosition = 5/12; // Center of third column (2/6 to 3/6)
    } else if (index === 3) {
      indexPosition = 7/12; // Center of fourth column (3/6 to 4/6)
    } else if (index === 4) {
      indexPosition = 9/12; // Center of fifth column (4/6 to 5/6)
    } else if (index === 5) {
      indexPosition = 11/12; // Center of sixth column (5/6 to 1)
    } else if (index === 6) {
      indexPosition = 1/12;
    } else if (index === 7) {
      indexPosition = 3/12;
    } else if (index === 8){
      indexPosition = 5/12;
    }
     else {
      indexPosition = 7/12;
    }
    return indexPosition;
  }

  const pieTitles0 = Object.entries(groupedHuman_sorted).map(
    ([key, element], index, array) => ({
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      yanchor: 'bottom',
      showarrow: false,
      index: index,
      text: (array[index][1][0].profile + array[index][1][0].matrix).padStart(
        8,
        ' '
      ),
      align: 'center',
      font: {
        weight: 'bold',
      },
      x: index < 4 ? 0.068 : index < 8 ? 0.206 : 0.344,
      y: index < 4 
        ? 1 - (index * 0.25) - 0.015
        : index < 8
        ? 1 - ((index - 4) * 0.25) - 0.015
        : 1 - ((index - 8) * 0.25) - 0.015,
    })
  );
  console.log("pieTitles0 ", pieTitles0);
  const pieTitles1 = Object.entries(groupedMm9_sorted).map(
    ([key, element], index, array) => ({
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      yanchor: 'bottom',
      showarrow: false,
      text: (array[index][1][0].profile + array[index][1][0].matrix).padStart(
        7,
        ' '
      ),
      align: 'center',
      font: {
        weight: 'bold',
      },
      x: 0.482, // Center of fourth column
      y: 1 - (index * 0.25) - 0.015,
    })
  );

  const pieTitles2 = Object.entries(groupedRn6).map(
    ([key, element], index, array) => ({
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      yanchor: 'bottom',
      showarrow: false,
      text: key.padStart(7, ' '),
      align: 'center',
      font: {
        weight: 'bold',
      },
      x: 0.620, // Center of fifth column
      y: 1 - (index * 0.25) - 0.015,
    })
  );

  const pieTitles3 = Object.entries(groupedGallus).map(
    ([key, element], index, array) => ({
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      yanchor: 'bottom',
      showarrow: false,
      text: key.padStart(7, ' '),
      align: 'center',
      font: {
        weight: 'bold',
      },
      x: 0.758, // Center of sixth column
      y: 1 - (index * 0.25) - 0.015,
    })
  );
  const pieTitles4 = Object.entries(groupCaenorhabditis).map(
    ([key, element], index, array) => ({
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      yanchor: 'bottom',
      showarrow: false,
      text: key.padStart(7, ' '),
      align: 'center',
      font: {
        weight: 'bold',
      },
      x: 0.896, // Center of seventh column
      y: 1 - (index * 0.25) - 0.01,
    })
  );
  const annotationTitle0 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'center',
    yanchor: 'bottom',
    showarrow: false,
    text: Object.keys(groupBySpecies)[0],
    font: {
      size: 14,
      weight: 'bold',
    },
    x: 0.206, // Center of all three columns
    y: 1.02,
  };

  const annotationTitle1 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'center',
    yanchor: 'bottom',
    showarrow: false,
    text: Object.keys(groupBySpecies)[1].replace(' ', '<br>'),
    font: {
      size: 12,
      weight: 'bold',
    },
    x: 0.482, // Center of fourth column
    y: 1.02,
  };

  const annotationTitle2 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'center',
    yanchor: 'bottom',
    showarrow: false,
    text: Object.keys(groupBySpecies)[2].replace(' ', '<br>'),
    font: {
      size: 12,
      weight: 'bold',
    },
    x: 0.620, // Center of fifth column
    y: 1.02,
  };

  const annotationTitle3 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'center',
    yanchor: 'bottom',
    showarrow: false,
    text: Object.keys(groupBySpecies)[3].replace(/\s+/g, '<br>'),
    font: {
      size: 12,
      weight: 'bold',
    },
    x: 0.758, // Center of sixth column
    y: 1.02,
  };

  const annotationTitle4 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'center',
    yanchor: 'bottom',
    showarrow: false,
    text: Object.keys(groupBySpecies)[4].replace(/\s+/g, '<br>'),
    font: {
      size: 12,
      weight: 'bold',
    },
    x: 0.896, // Center of seventh column
    y: 1.02,
  };

  const shapes = [
    // Vertical line between Homo sapiens and Mus musculus
    {
      type: 'line',
      xref: 'paper',
      yref: 'paper',
      x0: 0.412,
      y0: 0,
      x1: 0.412,
      y1: 1,
      line: {
        color: 'lightgray',
        width: 2,
        dash: 'solid',
      },
    },
    // Vertical line between Mus musculus and Rattus norvegicus
    {
      type: 'line',
      xref: 'paper',
      yref: 'paper',
      x0: 0.550,
      y0: 0,
      x1: 0.550,
      y1: 1,
      line: {
        color: 'lightgray',
        width: 2,
        dash: 'solid',
      },
    },
    // Vertical line between Rattus norvegicus and Gallus gallus
    {
      type: 'line',
      xref: 'paper',
      yref: 'paper',
      x0: 0.688,
      y0: 0,
      x1: 0.688,
      y1: 1,
      line: {
        color: 'lightgray',
        width: 2,
        dash: 'solid',
      },
    },
    // Vertical line between Gallus gallus and Caenorhabditis elegans
    {
      type: 'line',
      xref: 'paper',
      yref: 'paper',
      x0: 0.826,
      y0: 0,
      x1: 0.826,
      y1: 1,
      line: {
        color: 'lightgray',
        width: 2,
        dash: 'solid',
      },
    },
  ];

  const traces = [...tracePies0, ...tracePies1, ...tracePies2, ...tracePies3, ...tracePies4];
  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: 920,
    autosize: true,
    margin: {
      l: 0,
      r: 50,
      t: 100,
      b: 50,
    },
    legend: {
      title: {
        text: '\t <b>Signature Set Name</b>',
        font: {
          family: 'Times New Roman',
          size: 17,
        },
      },
      x: 1,
      xanchor: 'right',
      y: 0,
    },
    annotations: [
      ...pieTitles0,
      ...pieTitles1,
      ...pieTitles2,
      ...pieTitles3,
      ...pieTitles4,
      annotationTitle0,
      annotationTitle1,
      annotationTitle2,
      annotationTitle3,
      annotationTitle4
    ],
    shapes: [...shapes],
  };
  const config = {
    responsive: true,
    displayModeBar: true,
  };
  return { traces: traces, layout: layout, config };
}
