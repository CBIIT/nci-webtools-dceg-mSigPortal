export default function MSPrevalence(groupBySignature, groupBySample) {
  groupBySignature.sort((a, b) => b.samples.length - a.samples.length);
  console.log(groupBySignature);
  console.log(groupBySample);

  const tracesPie = {
    type: 'pie',
    labels: groupBySignature.map((group) => group.signatureName),
    values: groupBySignature.map((group) => group.samples.length),
    textposition: 'inside',
    textinfo: 'label+percent',
    showlegend: false,
  };
  console.log(tracesPie);

  const tracesBar = groupBySignature.map((group, groupIndex, array) => ({
    group: group,
    array: array,
    name: group.signatureName,
    type: 'bar',
    x: [group.signatureName],
    y: [group.samples.length / group.totalSamples],
    text: [
      Math.round((group.samples.length / group.totalSamples) * 1000) / 10 + '%',
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
  }));
  console.log(tracesBar);

  const traces = [tracesPie, ...tracesBar];

  console.log(traces);
  const layout = {
    grid: { rows: 1, columns: 2 },
    hoverlabel: { bgcolor: '#FFF' },
    height: 450,
    //width: 500,
    autosize: true,

    xaxis2: {
      showline: true,
      tickangle: -90,
      tickfont: {
        family: 'Arial, monospace',
      },
      tickmode: 'array',

      linecolor: '#E0E0E0',
      linewidth: 1,
      mirror: 'all',
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
      linecolor: '#D3D3D3',
      linewidth: 1,
      mirror: 'all',
      tickformat: ',.0%',
      showgrid: true,
      gridcolor: '#F5F5F5',
    },
  };

  return { traces: traces, layout: layout };
}
