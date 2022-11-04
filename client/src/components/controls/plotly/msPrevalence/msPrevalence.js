import { groupBy } from 'lodash';

export default function MSPrevalence(data, minimum) {
  // calculate median burden across cancer types
  const groupBySignature = groupBy(data, 'signatureName');
  console.log(data);
  console.log(groupBySignature);

  const dataResult = Object.entries(groupBySignature)
    .map(([signatureName, data]) => {
      const samples = data
        .filter((e) => e.exposure)
        .sort((a, b) => a.burden - b.burden);

      const burdens = samples.filter((e) => e.burden).map((e) => e.burden);

      const medianBurden =
        burdens.length % 2 == 0
          ? (burdens[burdens.length / 2] + burdens[burdens.length / 2 - 1]) / 2
          : burdens[Math.floor(burdens.length / 2)];

      return {
        signatureName,
        samples,
        medianBurden,
        totalSamples: data.length,
      };
    })
    .filter((e) => e.medianBurden)
    .sort((a, b) => a.medianBurden - b.medianBurden);
  // groupBySignature.sort(
  //   (a, b) =>
  //     a.samples.reduce((a, b) => a + b.exposure, 0) -
  //     b.samples.reduce((a, b) => a + b.exposure, 0)
  // );
  console.log(dataResult);
  let minumumNumber = 100;
  minimum === null || minimum === undefined
    ? (minumumNumber = 100)
    : (minumumNumber = parseInt(minimum));

  dataResult.sort(
    (a, b) =>
      b.samples.filter((e) => e.exposure >= minumumNumber).length /
        b.totalSamples -
      a.samples.filter((e) => e.exposure >= minumumNumber).length /
        a.totalSamples
  );
  const defaultNames = ['SBS', 'DBS', 'ID'];
  const names = dataResult.map((group) => group.signatureName);
  const longest = names.reduce((a, e) => (a > e.length ? a : e.length), 0);
  const extraMargin = longest < 10 ? 60 : longest * 7;
  console.log(names);
  console.log(longest);
  const contains = defaultNames.some((element) => {
    if (names[0].includes(element)) {
      return true;
    }

    return false;
  });
  console.log(contains);

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
  }

  const tracesPie = {
    type: 'pie',
    labels: dataResult.map((group) => group.signatureName),
    values: dataResult.map((group) =>
      group.samples.reduce((a, b) => a + b.exposure, 0)
    ),
    marker: {
      colors: dataResult.map((group) =>
        contains
          ? colors[group.signatureName.replace(/^\D*/, '').replace(')', '')]
          : colors[group.signatureName]
      ),
    },
    textposition: 'inside',
    texttemplate: '%{percent:.1%}',
    showlegend: false,
    marker: {
      colors: dataResult.map((group) =>
        contains
          ? colors[group.signatureName.replace(/^\D*/, '').replace(')', '')]
          : colors[group.signatureName]
      ),
    },
    test: dataResult.map((group) =>
      group.signatureName.replace(/^\D*/, '').replace(')', '')
    ),
    direction: 'clockwise',
    sort: true,
    domain: { x: [0, 0.2] },
    customdata: dataResult.map((e) => ({
      label: e.signatureName,
      value: e.samples.reduce((a, b) => a + b.exposure, 0),
    })),
    hovertemplate:
      '<b>Signature Name:</b>%{label}<br><b>Total sample: </b>%{value}<br>%{percent}<extra></extra>',
  };
  console.log(tracesPie);
  const tracesBar = dataResult.map((group, groupIndex, array) => ({
    name: group.signatureName,
    type: 'bar',
    marker: {
      color: contains
        ? colors[group.signatureName.replace(/^\D*/, '').replace(')', '')]
        : colors[group.signatureName],
    },
    x: [group.signatureName],
    y: [
      group.samples.filter((e) => e.exposure >= minumumNumber).length /
        group.totalSamples,
    ],
    text: [
      Math.round(
        (group.samples.filter((e) => e.exposure >= minumumNumber).length /
          group.totalSamples) *
          100 *
          10
      ) /
        10 +
        '%',
    ],
    customdata: [
      {
        signatureName: group.signatureName,
      },
    ],
    textposition: 'outside',
    xaxis: 'x2',
    yaxis: 'y2',
    //hoverinfo: 'x2+y2',
    hovertemplate:
      '<b> signatureName: ' +
      '</b>' +
      '%{customdata.signatureName}<br>' +
      '<b>Frequency: </b>' +
      '%{y:.1%}',
    showlegend: false,
    domain: {
      row: 0,
      column: 1,
    },
  }));

  console.log(tracesBar);
  const titleAnnotation = [
    {
      xref: 'paper',
      yref: 'paper',
      showarrow: false,
      x: 0.0225,
      y: 1.15,
      xanchor: 'top',
      text: '<b>Prevalence by mutations</b>',
      font: {
        size: 18,
        family: 'Arial',
      },
    },
    {
      xref: 'x2 domain',
      yref: 'paper',
      showarrow: false,
      x: 0.5,
      y: 1.15,
      xanchor: 'top',
      text: '<b>Prevalence by samples</b>',
      font: {
        size: 18,
        family: 'Arial',
      },
    },
  ];
  const barAnnotation = {
    xref: 'x2 domain',
    yref: 'paper',
    showarrow: false,
    x: 0.5,
    y: 0.5,
    xanchor: 'top',
    text: '<b>No signature with prevalence greater than 1%</b>',
    font: {
      size: 18,
      family: 'Arial',
    },
  };
  let titleAnnotations = [];

  const yMax = Math.max(...tracesBar.map((o) => o.y));
  let traces = [];
  if (yMax < 0.01) {
    traces = [tracesPie];
    titleAnnotations = [...titleAnnotation, barAnnotation];
  } else {
    traces = [tracesPie, ...tracesBar];
    titleAnnotations = [...titleAnnotation];
  }

  const layout = {
    grid: { rows: 1, columns: 2 },
    hoverlabel: { bgcolor: '#FFF' },
    height: 450,
    autosize: true,
    xaxis2: {
      showline: true,
      tickangle: -90,
      tickfont: {
        family: 'Arial, monospace',
      },
      tickmode: 'array',
      //tickvals: names.map((_, i) => i),
      //ticktext: names.map((e) => e.signatureName),
      linecolor: 'black',
      linewidth: 1,
      type: 'category',
      categoryorder: 'total descending',
      domain: [0.25, 1],
    },
    yaxis2: {
      title: {
        text: '<b>Frequency (%)</b>',
        font: {
          family: 'Times New Roman',
        },
      },

      range: [0, 1.1],
      ticks: 'inside',
      tickcolor: '#D3D3D3',
      linecolor: 'black',
      linewidth: 1,

      tickformat: ',.0%',
      showgrid: true,
      gridcolor: '#F5F5F5',
    },
    margin: {
      b: extraMargin,
    },
    annotations: titleAnnotations,
  };

  return { traces: traces, layout: layout };
}
