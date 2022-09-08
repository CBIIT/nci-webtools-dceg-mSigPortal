export default function MSPrevalence(groupBySignature, mutation) {
  groupBySignature.sort(
    (a, b) =>
      a.samples.reduce((a, b) => a + b.exposure, 0) -
      b.samples.reduce((a, b) => a + b.exposure, 0)
  );
  console.log(groupBySignature);
  let minumumNumber = 100;
  mutation === null || mutation === undefined
    ? (minumumNumber = 100)
    : (minumumNumber = parseInt(mutation));

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
    DBS1: '#4a9855',
    DBS2: '#e2a8ab',
    DBS3: '#40004b',
    DBS4: '#5aa1ca',
    DBS5: '#305d39',
    DBS6: '#785940',
    DBS7a: '#6e70b7',
    DBS7b: '#ff7f00',
    DBS7c: '#fec44f',
    DBS7d: '#846a2a',
    DBS8: '#cab2d6',
    DBS9: '#f4a582',
    DBS10a: '#8dd3c7',
    DBS10b: '#5e4fa2',
    DBS10c: '#761429',
    DBS11: '#9e0142',
    DBS12: '#ffed6f',
    DBS13: '#e41a1c',
    DBS14: '#ffffbf',
    DBS15: '#4d4d4d',
    DBS16: '#513276',
    DBS17a: '#df4c7d',
    DBS17b: '#08519c',
    DBS18: '#b3de69',
    DBS19: '#dfc27d',
    DBS20: '#b2182b',
    DBS21: '#9ecae1',
    DBS22: '#01665e',
    DBS23: '#d53e4f',
    DBS24: '#1c9099',
    DBS25: '#35978f',
    DBS26: '#ec7014',
    DBS27: '#f46d43',
    DBS28: '#de77ae',
    DBS29: '#fdae61',
    DBS30: '#d9d9d9',
    DBS31: '#f781bf',
    DBS32: '#dd1c77',
    DBS33: '#b25d7e',
    DBS34: '#fee08b',
    DBS35: '#fc8d59',
    DBS36: 'yellow',
    DBS37: '#e6f598',
    DBS38: '#abdda4',
    DBS39: '#636363',
    DBS40: '#b15928',
    DBS41: '#fccde5',
    DBS42: '#ae017e',
    DBS43: '#66c2a5',
    DBS44: '#8c6bb1',
    DBS45: '#3288bd',
    DBS46: '#e6f598',
    DBS47: '#bababa',
    DBS48: '#5e4fa2',
    DBS49: '#40004b',
    DBS50: '#762a83',
    DBS51: '#9970ab',
    DBS52: '#c2a5cf',
    DBS53: '#e7d4e8',
    DBS54: '#fcc5c0',
    DBS55: '#d9f0d3',
    DBS56: '#8c510a',
    DBS57: '#a6dba0',
    DBS58: '#5aae61',
    DBS59: '#1b7837',
    DBS60: '#00441b',
    DBS84: '#063C3C',
    DBS85: '#AA9139',
    DBS92: '#0E1844',
    'DBS-others': '#cececa',
    ID1: '#4a9855',
    ID2: '#e2a8ab',
    ID3: '#40004b',
    ID4: '#5aa1ca',
    ID5: '#305d39',
    ID6: '#785940',
    ID7a: '#6e70b7',
    ID7b: '#ff7f00',
    ID7c: '#fec44f',
    ID7d: '#846a2a',
    ID8: '#cab2d6',
    ID9: '#f4a582',
    ID10a: '#8dd3c7',
    ID10b: '#5e4fa2',
    ID10c: '#761429',
    ID11: '#9e0142',
    ID12: '#ffed6f',
    ID13: '#e41a1c',
    ID14: '#ffffbf',
    ID15: '#4d4d4d',
    ID16: '#513276',
    ID17a: '#df4c7d',
    ID17b: '#08519c',
    ID18: '#b3de69',
    ID19: '#dfc27d',
    ID20: '#b2182b',
    ID21: '#9ecae1',
    ID22: '#01665e',
    ID23: '#d53e4f',
    ID24: '#1c9099',
    ID25: '#35978f',
    ID26: '#ec7014',
    ID27: '#f46d43',
    ID28: '#de77ae',
    ID29: '#fdae61',
    ID30: '#d9d9d9',
    ID31: '#f781bf',
    ID32: '#dd1c77',
    ID33: '#b25d7e',
    ID34: '#fee08b',
    ID35: '#fc8d59',
    ID36: 'yellow',
    ID37: '#e6f598',
    ID38: '#abdda4',
    ID39: '#636363',
    ID40: '#b15928',
    ID41: '#fccde5',
    ID42: '#ae017e',
    ID43: '#66c2a5',
    ID44: '#8c6bb1',
    ID45: '#3288bd',
    ID46: '#e6f598',
    ID47: '#bababa',
    ID48: '#5e4fa2',
    ID49: '#40004b',
    ID50: '#762a83',
    ID51: '#9970ab',
    ID52: '#c2a5cf',
    ID53: '#e7d4e8',
    ID54: '#fcc5c0',
    ID55: '#d9f0d3',
    ID56: '#8c510a',
    ID57: '#a6dba0',
    ID58: '#5aae61',
    ID59: '#1b7837',
    ID60: '#00441b',
    ID84: '#063C3C',
    ID85: '#AA9139',
    ID92: '#0E1844',
    'ID-others': '#cececa',
  };

  const tracesPie = {
    type: 'pie',
    labels: groupBySignature.map((group) => group.signatureName),
    values: groupBySignature.map((group) =>
      group.samples.reduce((a, b) => a + b.exposure, 0)
    ),
    textposition: 'inside',
    texttemplate: '%{percent:.1%}',
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
    hovertemplate:
      '<b>Signature Name:</b>%{label}<br><b>Total sample: </b>%{value}<br>%{percent}<extra></extra>',
  };
  console.log(tracesPie);
  const tracesBar = groupBySignature.map((group, groupIndex, array) => ({
    group: group,
    array: array,
    name: group.signatureName,
    type: 'bar',
    minumumNumber: minumumNumber,
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
    customdata: groupBySignature.map((group) => ({
      signatureName: group.signatureName,
    })),
    textposition: 'outside',
    xaxis: 'x2',
    yaxis: 'y2',
    //hoverinfo: 'x2+y2',
    hovertemplate:
      '<b> signatureName: ' +
      '</b>' +
      '%{customdata.signatureName} <br>' +
      '<b>sample: </b>' +
      '%{y}',
    showlegend: false,
    domain: {
      row: 0,
      column: 1,
    },
    //hovertemplate: '<b>%{x}</b><br>%{y:.1%}<extra></extra>',
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
    text: '<b>No samples with prevalence greater than 1% </b>',
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

  const names = groupBySignature.map((group) => ({
    signatureName:
      group.signatureName.length > 10
        ? group.signatureName.substring(0, 10) + '...'
        : group.signatureName,
  }));
  console.log(names);
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
      tickvals: names.map((_, i) => i),
      ticktext: names.map((e) => e.signatureName),
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
