import { groupBy } from 'lodash';
export default function RsInMsigportal(rawData) {
  const profile_order_0 = [
    'ID29',
    'ID83',
    'CN48',
    'RS32',
    'RNA192',
    'SBS96',
    'SBS192',
    'SBS288',
    'SBS1536',
    'DBS78',
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

  const groupedHuman_sorted = profile_order_0.map((key) => groupedHuman[key]).slice(0, 5);
  console.log("groupedHuman_sorted ", groupedHuman_sorted);
  const groupedMm9_sorted = profile_order_1.map((key) => groupedMm9[key]);
  
  // Create additional groups for the remaining Human data (5 more for column 2)
  const groupedHuman_sorted_col2 = profile_order_0.map((key) => groupedHuman[key]).slice(5, 10);
  
  const tracePies0 = Object.entries(groupedHuman_sorted).map(
    ([key, element], index, array) => ({
      type: 'pie',
      marker: {
        color: element ? element.map((e) => colors[e.signatureSetName]) : [],
        line: {
          color: 'black',
          width: 1,
        },
      },
      textposition: 'inside',
      labels: element ? element.map((e) => e.signatureSetName) : [],
      values: element ? element.map((e) => parseInt(e.count)) : [],
      texttemplate: '%{value}',
      direction: 'clockwise',
      index: index,
      name: element && element[0] ? element[0].profile + element[0].matrix : '',
      domain: {
        x: [0, 0.15],
        y: [
          Math.round((0.84 - index * 0.18) * 100) / 100,
          Math.round((0.84 - index * 0.18 + 0.16) * 100) / 100,
        ],
      },
      hovertemplate:
        '<b>%{label}</b> <br>%{percent} </br> %{value}  <extra></extra>',
    })
  ).filter(pie => pie.name !== ''); // Filter out empty pies
  console.log("tracePies0 ", tracePies0);
  
  // Column 2: Remaining 5 Human pies
  const tracePies0_col2 = Object.entries(groupedHuman_sorted_col2).map(
    ([key, element], index, array) => ({
      type: 'pie',
      marker: {
        color: element ? element.map((e) => colors[e.signatureSetName]) : [],
        line: {
          color: 'black',
          width: 1,
        },
      },
      textposition: 'inside',
      labels: element ? element.map((e) => e.signatureSetName) : [],
      values: element ? element.map((e) => parseInt(e.count)) : [],
      texttemplate: '%{value}',
      direction: 'clockwise',
      index: index,
      name: element && element[0] ? element[0].profile + element[0].matrix : '',
      domain: {
        x: [0.17, 0.32],
        y: [
          Math.round((0.84 - index * 0.18) * 100) / 100,
          Math.round((0.84 - index * 0.18 + 0.16) * 100) / 100,
        ],
      },
      hovertemplate:
        '<b>%{label}</b> <br>%{percent} </br> %{value}  <extra></extra>',
    })
  ).filter(pie => pie.name !== '');
  
  // Column 3: Mouse data (3 pies)
  const tracePies1 = Object.entries(groupedMm9_sorted).map(
    ([key, element], index, array) => ({
      type: 'pie',
      e: element,
      marker: {
        color: element ? element.map((e) => colors[e.signatureSetName]) : [],
        line: {
          color: 'black',
          width: 1,
        },
      },
      textposition: 'inside',
      labels: element ? element.map((e) => e.signatureSetName) : [],
      values: element ? element.map((e) => parseInt(e.count)) : [],
      texttemplate: '%{value}',
      direction: 'clockwise',
      name: key,
      domain: {
        x: [0.34, 0.49],
        y: [
          Math.round((0.84 - index * 0.18) * 100) / 100,
          Math.round((0.84 - index * 0.18 + 0.16) * 100) / 100,
        ],
      },
      hovertemplate:
        '<b>%{label}</b> <br>%{percent} </br> %{value}  <extra></extra>',
    })
  ).filter(pie => pie.labels && pie.labels.length > 0);

  // Column 4: Rat data (2 pies)
  const tracePies2 = Object.entries(groupedRn6).slice(0, 2).map(
    ([key, element], index, array) => ({
      type: 'pie',
      e: element,
      marker: {
        color: element ? element.map((e) => colors[e.signatureSetName]) : [],
        line: {
          color: 'black',
          width: 1,
        },
      },
      textposition: 'inside',
      labels: element ? element.map((e) => e.signatureSetName) : [],
      values: element ? element.map((e) => parseInt(e.count)) : [],
      texttemplate: '%{value}',
      direction: 'clockwise',
      name: key,
      domain: {
        x: [0.51, 0.66],
        y: [
          Math.round((0.84 - index * 0.18) * 100) / 100,
          Math.round((0.84 - index * 0.18 + 0.16) * 100) / 100,
        ],
      },
      hovertemplate:
        '<b>%{label}</b> <br>%{percent} </br> %{value}  <extra></extra>',
    })
  ).filter(pie => pie.labels && pie.labels.length > 0);
  // Column 5: Gallus data (1 pie)
  const tracePies3 = Object.entries(groupedGallus).slice(0, 1).map(
    ([key, element], index, array) => ({
      type: 'pie',
      e: element,
      marker: {
        color: element ? element.map((e) => colors[e.signatureSetName]) : [],
        line: {
          color: 'black',
          width: 1,
        },
      },
      textposition: 'inside',
      labels: element ? element.map((e) => e.signatureSetName) : [],
      values: element ? element.map((e) => parseInt(e.count)) : [],
      texttemplate: '%{value}',
      direction: 'clockwise',
      name: key,
      domain: {
        x: [0.68, 0.83],
        y: [0.74, 0.90],
      },
      hovertemplate:
        '<b>%{label}</b> <br>%{percent} </br> %{value}  <extra></extra>',
    })
  ).filter(pie => pie.labels && pie.labels.length > 0);

  // Column 6: Caenorhabditis data (1 pie)
  const tracePies4 = Object.entries(groupCaenorhabditis).slice(0, 1).map(
    ([key, element], index, array) => ({
      type: 'pie',
      e: element,
      marker: {
        color: element ? element.map((e) => colors[e.signatureSetName]) : [],
        line: {
          color: 'black',
          width: 1,
        },
      },
      textposition: 'inside',
      labels: element ? element.map((e) => e.signatureSetName) : [],
      values: element ? element.map((e) => parseInt(e.count)) : [],
      texttemplate: '%{value}',
      direction: 'clockwise',
      name: key,
      domain: {
        x: [0.85, 1.0],
        y: [0.74, 0.90],
      },
      hovertemplate:
        '<b>%{label}</b> <br>%{percent} </br> %{value}  <extra></extra>',
    })
  ).filter(pie => pie.labels && pie.labels.length > 0);

  function indexPos(index, column = 0) {
    const columnPositions = [0.09, 0.29, 0.49, 0.69, 0.89]; // 5 columns for Human
    const column2Positions = [0.29, 0.29, 0.29]; // Mouse column
    const singleColumnPosition = 0.49; // For other species with fewer charts
    
    if (column === 0) { // Human - 5 columns
      return columnPositions[index % 5];
    } else if (column === 1) { // Mouse - single column
      return 0.29;
    } else { // Other species - single column each
      if (column === 2) return 0.49; // Rat
      if (column === 3) return 0.69; // Gallus
      if (column === 4) return 0.89; // Caenorhabditis
    }
    return 0.5;
  }

  // Column 1 titles
  const pieTitles0 = Object.entries(groupedHuman_sorted).map(
    ([key, element], index, array) => ({
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      yanchor: 'middle',
      showarrow: false,
      index: index,
      text: element && element[0] ? (element[0].profile + element[0].matrix).padStart(8, ' ') : '',
      align: 'center',
      x: 0.075,
      y: Math.round((0.84 - index * 0.18 + 0.08) * 100) / 100,
    })
  ).filter(title => title.text !== '');

  // Column 2 titles  
  const pieTitles0_col2 = Object.entries(groupedHuman_sorted_col2).map(
    ([key, element], index, array) => ({
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      yanchor: 'middle',
      showarrow: false,
      index: index,
      text: element && element[0] ? (element[0].profile + element[0].matrix).padStart(8, ' ') : '',
      align: 'center',
      x: 0.245,
      y: Math.round((0.84 - index * 0.18 + 0.08) * 100) / 100,
    })
  ).filter(title => title.text !== '');
  
  console.log("pieTitles0 ", pieTitles0);
  
  // Column 3 titles
  const pieTitles1 = Object.entries(groupedMm9_sorted).map(
    ([key, element], index, array) => ({
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      yanchor: 'middle',
      showarrow: false,
      text: element && element[0] ? (element[0].profile + element[0].matrix).padStart(7, ' ') : '',
      align: 'center',
      x: 0.415,
      y: Math.round((0.84 - index * 0.18 + 0.08) * 100) / 100,
    })
  ).filter(title => title.text !== '');

  // Column 4 titles
  const pieTitles2 = Object.entries(groupedRn6).slice(0, 2).map(
    ([key, element], index, array) => ({
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      yanchor: 'middle',
      showarrow: false,
      text: key.padStart(7, ' '),
      align: 'center',
      x: 0.585,
      y: Math.round((0.84 - index * 0.18 + 0.08) * 100) / 100,
    })
  );

  // Column 5 titles
  const pieTitles3 = Object.entries(groupedGallus).slice(0, 1).map(
    ([key, element], index, array) => ({
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      yanchor: 'top',
      showarrow: false,
      text: key.padStart(7, ' '),
      align: 'center',
      x: 0.755,
      y: 0.82,
    })
  );
  
  // Column 6 titles
  const pieTitles4 = Object.entries(groupCaenorhabditis).slice(0, 1).map(
    ([key, element], index, array) => ({
      xref: 'paper',
      yref: 'paper',
      xanchor: 'center',
      yanchor: 'top',
      showarrow: false,
      text: key.padStart(7, ' '),
      align: 'center',
      x: 0.925,
      y: 0.82,
    })
  );
  const annotationTitle0 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'center',
    yanchor: 'bottom',
    showarrow: false,
    text: Object.keys(groupBySpecies)[0] + ' (1-5)',
    font: {
      size: 14,
    },
    x: 0.075,
    y: 1.02,
  };

  const annotationTitle1 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'center',
    yanchor: 'bottom',
    showarrow: false,
    text: Object.keys(groupBySpecies)[0] + ' (6-10)',
    font: {
      size: 14,
    },
    x: 0.245,
    y: 1.02,
  };

  const annotationTitle2 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'center',
    yanchor: 'bottom',
    showarrow: false,
    text: Object.keys(groupBySpecies)[1],
    font: {
      size: 14,
    },
    x: 0.415,
    y: 1.02,
  };

  const annotationTitle3 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'center',
    yanchor: 'bottom',
    showarrow: false,
    text: Object.keys(groupBySpecies)[2],
    font: {
      size: 14,
    },
    x: 0.585,
    y: 1.02,
  };

  const annotationTitle4 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'center',
    yanchor: 'top',
    showarrow: false,
    text: Object.keys(groupBySpecies)[3],
    font: {
      size: 14,
    },
    x: 0.755,
    y: 0.72,
  };

  const annotationTitle5 = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'center',
    yanchor: 'top',
    showarrow: false,
    text: Object.keys(groupBySpecies)[4],
    font: {
      size: 14,
    },
    x: 0.925,
    y: 0.72,
  };

  const shapes = [
    // Vertical lines separating columns
    {
      type: 'line',
      xref: 'paper',
      yref: 'paper',
      x0: 0.16,
      y0: 0,
      x1: 0.16,
      y1: 1.0,
      line: {
        color: 'lightgray',
        width: 1,
      },
    },
    {
      type: 'line',
      xref: 'paper',
      yref: 'paper',
      x0: 0.33,
      y0: 0,
      x1: 0.33,
      y1: 1.0,
      line: {
        color: 'lightgray',
        width: 1,
      },
    },
    {
      type: 'line',
      xref: 'paper',
      yref: 'paper',
      x0: 0.50,
      y0: 0,
      x1: 0.50,
      y1: 1.0,
      line: {
        color: 'lightgray',
        width: 1,
      },
    },
    {
      type: 'line',
      xref: 'paper',
      yref: 'paper',
      x0: 0.67,
      y0: 0,
      x1: 0.67,
      y1: 1.0,
      line: {
        color: 'lightgray',
        width: 1,
      },
    },
    {
      type: 'line',
      xref: 'paper',
      yref: 'paper',
      x0: 0.84,
      y0: 0,
      x1: 0.84,
      y1: 1.0,
      line: {
        color: 'lightgray',
        width: 1,
      },
    },
  ];

  const traces = [...tracePies0, ...tracePies0_col2, ...tracePies1, ...tracePies2, ...tracePies3, ...tracePies4];
  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: 800,
    width: 1600,
    autosize: true,
    legend: {
      title: {
        text: '\t <b>Signature Set Name</b>',
        font: {
          family: 'Times New Roman',
          size: 17,
        },
      },
      x: 0.84,
      xanchor: 'left',
      y: 0.6,
      yanchor: 'top',
    },
    annotations: [
      ...pieTitles0,
      ...pieTitles0_col2,
      ...pieTitles1,
      ...pieTitles2,
      ...pieTitles3,
      ...pieTitles4,
      annotationTitle0,
      annotationTitle1,
      annotationTitle2,
      annotationTitle3,
      annotationTitle4,
      annotationTitle5
    ],
    shapes: [...shapes],
  };
  const config = {
    //responsive: true,
  };
  return { traces: traces, layout: layout, config };
}
