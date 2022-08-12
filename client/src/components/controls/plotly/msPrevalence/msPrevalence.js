export default function MSPrevalence(groupBySignature, mutation) {
  groupBySignature.sort(
    (a, b) =>
      a.samples.reduce((a, b) => a + b.exposure, 0) -
      b.samples.reduce((a, b) => a + b.exposure, 0)
  );

  let minumumNumber = 100;
  mutation === null || mutation === undefined
    ? (minumumNumber = 100)
    : (minumumNumber = parseInt(mutation));

  // const groupBySignature_sortExposure = groupBySignature.sort(
  //   (a, b) =>
  //     b.samples.filter((e) => e.exposure >= minumumNumber).length /
  //       a.totalSamples -
  //     a.samples.filter((e) => e.exposure >= minumumNumber).length /
  //       b.totalSamples
  // );

  // console.log(groupBySignature_sortExposure);
  const colors = {
    SBS1: '#4a9855',
    SBS2: '#e2a8ab',
    SBS3: '#40004b',
    SBS4: '#5aa1ca',
    SBS5: '#305d39',
    SBS6: '#785940',
    SBS7a: '#6e70b7',
    SBS7b: '#ff7f00',
    SBS7c: '#fec44f',
    SBS7d: '#846a2a',
    SBS8: '#cab2d6',
    SBS9: '#f4a582',
    SBS10a: '#8dd3c7',
    SBS10b: '#5e4fa2',
    SBS10c: '#761429',
    SBS11: '#9e0142',
    SBS12: '#ffed6f',
    SBS13: '#e41a1c',
    SBS14: '#ffffbf',
    SBS15: '#4d4d4d',
    SBS16: '#513276',
    SBS17a: '#df4c7d',
    SBS17b: '#08519c',
    SBS18: '#b3de69',
    SBS19: '#dfc27d',
    SBS20: '#b2182b',
    SBS21: '#9ecae1',
    SBS22: '#01665e',
    SBS23: '#d53e4f',
    SBS24: '#1c9099',
    SBS25: '#35978f',
    SBS26: '#ec7014',
    SBS27: '#f46d43',
    SBS28: '#de77ae',
    SBS29: '#fdae61',
    SBS30: '#d9d9d9',
    SBS31: '#f781bf',
    SBS32: '#dd1c77',
    SBS33: '#b25d7e',
    SBS34: '#fee08b',
    SBS35: '#fc8d59',
    SBS36: 'yellow',
    SBS37: '#e6f598',
    SBS38: '#abdda4',
    SBS39: '#636363',
    SBS40: '#b15928',
    SBS41: '#fccde5',
    SBS42: '#ae017e',
    SBS43: '#66c2a5',
    SBS44: '#8c6bb1',
    SBS45: '#3288bd',
    SBS46: '#e6f598',
    SBS47: '#bababa',
    SBS48: '#5e4fa2',
    SBS49: '#40004b',
    SBS50: '#762a83',
    SBS51: '#9970ab',
    SBS52: '#c2a5cf',
    SBS53: '#e7d4e8',
    SBS54: '#fcc5c0',
    SBS55: '#d9f0d3',
    SBS56: '#8c510a',
    SBS57: '#a6dba0',
    SBS58: '#5aae61',
    SBS59: '#1b7837',
    SBS60: '#00441b',
    SBS84: '#063C3C',
    SBS85: '#AA9139',
    SBS92: '#0E1844',
    'SBS-others': '#cececa',
  };

  const tracesPie = {
    type: 'pie',
    labels: groupBySignature.map((group) => group.signatureName),
    values: groupBySignature.map((group) =>
      group.samples.reduce((a, b) => a + b.exposure, 0)
    ),
    textposition: 'inside',
    textinfo: 'percent',
    showlegend: false,
    marker: {
      colors: groupBySignature.map((group) => colors[group.signatureName]),
    },
    direction: 'clockwise',
    sort: false,
    domain: { x: [0, 0.2] },
    customdata: groupBySignature.map((e) => ({
      label: e.signatureName,
      value: e.samples.reduce((a, b) => a + b.exposure, 0),
    })),
    hovertemplate: '<b>%{label}</b><br>%{value}<br>%{percent}<extra></extra>',
  };
  console.log(tracesPie);
  const tracesBar = groupBySignature.map((group, groupIndex, array) => ({
    group: group,
    array: array,
    name: group.signatureName,
    type: 'bar',
    minumumNumber: minumumNumber,
    lenght: group.samples.filter((e) => e.exposure >= 100),
    marker: {
      color: colors[group.signatureName],
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
    textposition: 'outside',
    xaxis: 'x2',
    yaxis: 'y2',
    hoverinfo: 'x2+y2',
    showlegend: false,
    domain: {
      row: 0,
      column: 1,
    },
    hovertemplate: '<b>%{x}</b><br>%{y:.1%}<extra></extra>',
  }));

  const traces = [tracesPie, ...tracesBar];

  const titleAnnotations = [
    {
      xref: 'paper',
      yref: 'paper',
      showarrow: false,
      x: 0.05,
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

  const layout = {
    grid: { rows: 1, columns: 2 },
    hoverlabel: { bgcolor: '#FFF' },
    height: 450,

    autosize: true,
    title: '<b>Prevalence of Mutational Signatures</b>',

    xaxis2: {
      showline: true,
      tickangle: -90,
      tickfont: {
        family: 'Arial, monospace',
      },
      tickmode: 'array',
      linecolor: 'black',
      linewidth: 1,
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
    annotations: titleAnnotations,
  };

  return { traces: traces, layout: layout };
}
